const VSHADER_SOURCE = `
  attribute vec4 a_Position;
  attribute vec4 a_Color;
  attribute vec4 a_Normal;
  uniform mat4 u_MvpMatrix;
  uniform mat4 u_ModelMatrix; // Model matrix
  uniform mat4 u_NormalMatrix; // Transformation matrix of the normal
  varying vec4 v_Color;
  varying vec3 v_Normal;
  varying vec3 v_Position;
  void main() {
    gl_PointSize = 3.5;
    vec4 color = a_Color;
    gl_Position = u_MvpMatrix * a_Position;
    // Calculate the vertex position in the world coordinate
    v_Position = vec3(u_ModelMatrix * a_Position);
    v_Normal = normalize(vec3(u_NormalMatrix * a_Normal));
    v_Color = color;
  }`;

const FSHADER_SOURCE = `
  #ifdef GL_ES
  precision mediump float;
  #endif
  uniform vec3 u_LightColor; // Light color
  uniform vec3 u_LightPosition; // Position of the light source
  uniform vec3 u_AmbientLight; // Ambient light color
  varying vec3 v_Normal;
  varying vec3 v_Position;
  varying vec4 v_Color;
  void main() {
    // Normalize the normal because it is interpolated and not 1.0 in length any more
    vec3 normal = normalize(v_Normal);
    // Calculate the light direction and make it 1.0 in length
    vec3 lightDirection = normalize(u_LightPosition - v_Position);
    // The dot product of the light direction and the normal
    float nDotL = max(dot(lightDirection, normal), 0.0);
    // Calculate the final color from diffuse reflection and ambient reflection
    vec3 diffuse = u_LightColor * v_Color.rgb * nDotL;
    vec3 ambient = u_AmbientLight * v_Color.rgb;
    gl_FragColor = vec4(diffuse + ambient, v_Color.a);
  }`;


class Camera {
  constructor() {
    this.pos = DEFAULT_CAMERA_POS;
    this.lookAt = DEFAULT_CAMERA_LOOK_AT;
    this.up = DEFAULT_CAMERA_UP;
  }
}


class Graphic {
  constructor() {
    this.canvas = document.getElementById('webgl');

    this.camera = new Camera();
    this.fpsInterval = 1000 / FPS_LIMIT;
    this.lastUpdateTime = Date.now();

    /* params needed for scene rotation by mouse movement */
    this.drag = false;
    this.previousX;
    this.previousY;
    this.dX = 0;
    this.dY = 0;
    this.THETA = -10;
    this.PHI = 6;

    this.canvas.onmousedown = e => {
      this.drag = true;
      this.previousX = e.pageX;
      this.previousY = e.pageY;
      e.preventDefault();
      return false;
    };
    this.canvas.onmouseup = _ => {
      this.drag = false;
    };
    this.canvas.onmousemove = e => {
      if (!this.drag) return false;
      this.dX = (e.pageX - this.previousX) * 2 * Math.PI / 10;
      this.dY = (e.pageY - this.previousY) * 2 * Math.PI / 10;

      this.THETA += this.dX;
      this.PHI += this.dY;

      this.THETA %= 360;
      this.PHI %= 360;

      this.previousX = e.pageX;
      this.previousY = e.pageY;
      e.preventDefault();
    };
    /* end of params needed for scene rotation */

    this.gl = getWebGLContext(this.canvas);
    if (!this.gl) {
      console.error('Failed to get the rendering context for WebGL');
      return;
    }

    if (!initShaders(this.gl, VSHADER_SOURCE, FSHADER_SOURCE)) {
      console.error('Failed to intialize shaders.');
      return;
    }

    this.u_ModelMatrix = this.gl.getUniformLocation(this.gl.program, 'u_ModelMatrix');
    this.u_MvpMatrix = this.gl.getUniformLocation(this.gl.program, 'u_MvpMatrix');
    this.u_NormalMatrix = this.gl.getUniformLocation(this.gl.program, 'u_NormalMatrix');

    this.u_LightColor = this.gl.getUniformLocation(this.gl.program, 'u_LightColor');
    this.u_LightPosition = this.gl.getUniformLocation(this.gl.program, 'u_LightPosition');
    this.u_AmbientLight = this.gl.getUniformLocation(this.gl.program, 'u_AmbientLight');

    if (!this.u_ModelMatrix || !this.u_NormalMatrix || !this.u_MvpMatrix || !this.u_LightColor || !this.u_LightPosition || !this.u_AmbientLight) {
      console.error('Failed to get the storage location');
      return;
    }

    this.gl.uniform3f(this.u_LightColor, ...LIGHT_COLOR);
    this.gl.uniform3f(this.u_LightPosition, ...LIGHT_POS);
    this.gl.uniform3f(this.u_AmbientLight, ...AMBIENT_LIGHT);

    this.gl.enable(this.gl.DEPTH_TEST);
  }

  setupSphere(sphere) {
    this.sphere = sphere;
  }

  setupPoints(positions) {
    this.points = {
      positions: positions,
      indices: [...Array(positions.length / 3).keys()],
      color: COLORS.BLACK,
      mode: DRAW_MODE.POINTS
    };
  }

  setupIsolines(isolines) {
    this.isolines = [];

    for (let isoline of isolines) {
      if (isoline.length > MAX_POINTS_IN_ISOLINE) {
        // simple technique for better performance
        let k = Math.ceil(isoline.length / MAX_POINTS_IN_ISOLINE);
        isoline = isoline.filter((_, i) => i % k === 0);
      }

      isoline = isoline.flat();

      this.isolines.push({
        positions: isoline,
        indices: [...Array(isoline.length / 3).keys()],
        color: COLORS.BLUE,
        mode: DRAW_MODE.LINE_LOOP
      });
    }
  }

  setupDebugIsolines(isolines) {
    this.debugIsolines = [];

    for (let isoline of isolines) {
      if (isoline.length > MAX_POINTS_IN_ISOLINE) {
        // simple technique for better performance
        let k = Math.ceil(isoline.length / MAX_POINTS_IN_ISOLINE);
        isoline = isoline.filter((_, i) => i % k === 0);
      }

      isoline = isoline.flat();

      this.debugIsolines.push({
        positions: isoline,
        indices: [...Array(isoline.length / 3).keys()],
        color: COLORS.RED,
        mode: DRAW_MODE.LINE_LOOP
      });
    }
  }

  setupMeanPoint(point) {
    let pointLen = Math.hypot(...point);
    let normedPoint = point.map(x => x / pointLen);
    this.meanPoint = {
      positions: normedPoint,
      indices: [0],
      color: COLORS.BLUE,
      mode: DRAW_MODE.POINTS
    };
  }

  setupCoordLines() {
    let delta = 0.015;
    let positions = [
      // X coordinate line with an arrow
      1.3,0,0, 1.24,delta,delta, 1.24,-delta,delta, 1.3,0,0,
      1.24,delta,delta, 1.24,delta,-delta, 1.3,0,0,
      1.24,-delta,-delta, 1.24,-delta,delta, 1.3,0,0,
      1.24,-delta,-delta, 1.24,delta,-delta, 1.3,0,0,
      0,0,0,
      // Y coordinate line with an arrow
      0,1.3,0, delta,1.24,delta, -delta,1.24,delta, 0,1.3,0,
      -delta,1.24,delta, -delta,1.24,-delta, 0,1.3,0,
      -delta,1.24,-delta, delta,1.24,-delta, 0,1.3,0,
      delta,1.24,-delta, delta,1.24,delta, 0,1.3,0,
      0,0,0,
      // Z coordinate line with an arrow
      0,0,1.3, delta,delta,1.24, delta,-delta,1.24, 0,0,1.3,
      delta,-delta,1.24, -delta,-delta,1.24, 0,0,1.3,
      -delta,-delta,1.24, -delta,delta,1.24, 0,0,1.3,
      -delta,delta,1.24, delta,delta,1.24, 0,0,1.3,
    ];

    this.coordLines = {
      positions: positions,
      indices: [...Array(positions.length / 3).keys()],
      color: COLORS.BLUE,
      mode: DRAW_MODE.LINE_STRIP
    };
  }

  // these points show the start from which s-statistic is counted
  setupDebugIntPoints(positions) {
    this.debugIntPoints = {
      positions: positions,
      indices: [...Array(positions.length / 3).keys()],
      color: COLORS.BLUE,
      mode: DRAW_MODE.POINTS
    };
  }

  // these points show in which direction from debugIntPoints s-stat is counted
  setupDebugDirPoints(positions) {
    this.debugDirPoints = {
      positions: positions,
      indices: [...Array(positions.length / 3).keys()],
      color: COLORS.LIGHT_BLUE,
      mode: DRAW_MODE.POINTS
    };
  }

  initVertexBuffersForObject(object) {
    if (!this.initArrayBuffer('a_Position', new Float32Array(object.positions), this.gl.FLOAT, 3)) return false;
    if (!this.initArrayBuffer('a_Normal', new Float32Array(object.positions), this.gl.FLOAT, 3)) return false;

    let colorArray = new Float32Array(object.indices.length * 4);
    for (let i = 0; i < colorArray.length; ++i) {
      colorArray[i] = object.color[i % 4];
    }
    if (!this.initArrayBuffer('a_Color', colorArray, this.gl.FLOAT, 4)) return false;

    let indexBuffer = this.gl.createBuffer();
    if (!indexBuffer) {
      console.error('Failed to create the buffer object');
      return false;
    }
    this.gl.bindBuffer(this.gl.ELEMENT_ARRAY_BUFFER, indexBuffer);
    this.gl.bufferData(this.gl.ELEMENT_ARRAY_BUFFER, new Uint32Array(object.indices), this.gl.STATIC_DRAW);

    return true;
  }

  initArrayBuffer(attribute, data, type, num) {
    let buffer = this.gl.createBuffer();
    if (!buffer) {
      console.error('Failed to create the buffer object');
      return false;
    }

    this.gl.bindBuffer(this.gl.ARRAY_BUFFER, buffer);
    this.gl.bufferData(this.gl.ARRAY_BUFFER, data, this.gl.STATIC_DRAW);

    let a_attribute = this.gl.getAttribLocation(this.gl.program, attribute);
    if (a_attribute < 0) {
      console.error('Failed to get the storage location of ' + attribute);
      return false;
    }

    this.gl.vertexAttribPointer(a_attribute, num, type, false, 0, 0);
    this.gl.enableVertexAttribArray(a_attribute);
    this.gl.bindBuffer(this.gl.ARRAY_BUFFER, null);

    return true;
  }

  clear() {
    this.gl.clearColor(1, 1, 1, 1)
    this.gl.clear(this.gl.COLOR_BUFFER_BIT | this.gl.DEPTH_BUFFER_BIT);
  }

  // object must have positions array with 3-dimensional coordinates, indices array, color array and mode property (1 of DRAW_MODE constant)
  draw(object) {
    this.gl.clear(this.gl.DEPTH_BUFFER_BIT)

    if (!this.initVertexBuffersForObject(object)) {
      console.error('Failed to set the vertex information');
      return;
    }

    let modelMatrix = new Matrix4();
    let mvpMatrix = new Matrix4();
    let normalMatrix = new Matrix4();

    modelMatrix.setRotate(this.THETA, 0, 0, 1); // Rotate around the y-axis
    modelMatrix.rotate(this.PHI, 0, 1, 0); // Rotate around the x-axis

    mvpMatrix.setPerspective(30, this.canvas.width / this.canvas.height, 1, 100);
    mvpMatrix.lookAt(
      this.camera.pos.x, this.camera.pos.y, this.camera.pos.z,
      this.camera.lookAt.x, this.camera.lookAt.y, this.camera.lookAt.z,
      this.camera.up.x, this.camera.up.y, this.camera.up.z);
    mvpMatrix.multiply(modelMatrix);

    normalMatrix.setInverseOf(modelMatrix);
    normalMatrix.transpose();

    this.gl.uniformMatrix4fv(this.u_ModelMatrix, false, modelMatrix.elements);
    this.gl.uniformMatrix4fv(this.u_MvpMatrix, false, mvpMatrix.elements);
    this.gl.uniformMatrix4fv(this.u_NormalMatrix, false, normalMatrix.elements);

    this.gl.getExtension('OES_element_index_uint'); // for UNSIGNED_INT support

    let vertex_n = object.indices.length
    switch(object.mode) {
    case DRAW_MODE.POINTS:
      this.gl.drawElements(this.gl.POINTS, vertex_n, this.gl.UNSIGNED_INT, 0);
      break;
    case DRAW_MODE.LINE_STRIP:
      this.gl.drawElements(this.gl.LINE_STRIP, vertex_n, this.gl.UNSIGNED_INT, 0);
      break;
    case DRAW_MODE.LINE_LOOP:
      this.gl.drawElements(this.gl.LINE_LOOP, vertex_n, this.gl.UNSIGNED_INT, 0);
      break;
    case DRAW_MODE.TRIANGLE_STRIP:
      this.gl.drawElements(this.gl.TRIANGLE_STRIP, vertex_n, this.gl.UNSIGNED_INT, 0);
      break;
    }
  }

  drawScene() {
    this.clear();

    if (this.coordLines) {
      this.draw(this.coordLines);
    }
    if (this.sphere) {
      this.draw(this.sphere);
    }
    if (this.isolines && this.isolines.length) {
      this.isolines.forEach(line => this.draw(line));
    }
    if (this.debugIsolines && this.debugIsolines.length) {
      this.debugIsolines.forEach(line => this.draw(line));
    }
    if (this.meanPoint) {
      this.draw(this.meanPoint);
    }
    if (this.debugIntPoints) { // TODO (debug)
      this.draw(this.debugIntPoints);
    }
    if (this.debugDirPoints) { // TODO (debug)
      this.draw(this.debugDirPoints);
    }
    if (this.points) {
      this.draw(this.points);
    }
  }

  tick() {
    requestAnimationFrame(_ => this.tick());

    let time = Date.now();
    let elapsed = time - this.lastUpdateTime;
    if (elapsed > this.fpsInterval) {
      this.drawScene();
      this.lastUpdateTime = time - (elapsed % this.fpsInterval);
    }
  }
}
