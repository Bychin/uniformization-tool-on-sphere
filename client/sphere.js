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
    gl_PointSize = 3.5; // TODO
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


const SPIRAL_SPHERE = 'spiral';
const CLASSIC_SPHERE = 'classic';

const SPHERE_COLOR = [1., 0.2, 0.2, 0.5];
const POINTS_COLOR = [0., 0., 0., 1.];

const LIGHT_COLOR = [0.8, 0.8, 0.8];
const LIGHT_POS = [-5.0, 0.0, -2.0];
const AMBIENT_LIGHT = [0.2, 0.2, 0.2];

const DEFAULT_SPHERE_ANGLE_DEG = 0;


const DRAW_MODE = {
  POINTS: 'points',
  LINE_STRIP: 'line_strip',
  TRIANGLE_STRIP: 'triangle_strip'
}


const DEFAULT_CAMERA_POS = {
  x: 5,
  y: 0,
  z: 0
}
const DEFAULT_CAMERA_LOOK_AT = {
  x: 0,
  y: 0,
  z: 0
}
const DEFAULT_CAMERA_UP = {
  x: 0,
  y: 0,
  z: 1
}

const AMORTIZATION = 0.;


class Sphere {
  // TODO others angles of sphere
  // TODO handler on mouse drag

  constructor(type, mode, div, radius, center) {
    this.mode = mode;
    this.color = SPHERE_COLOR;

    this.DIV = div;
    this.R = radius;
    this.center = center;
    this.angleZ = DEFAULT_SPHERE_ANGLE_DEG;

    if (type === SPIRAL_SPHERE) {
      this.type = SPIRAL_SPHERE;
      this.generateSpiralSphere();

    } else if (type === CLASSIC_SPHERE) {
      this.type = CLASSIC_SPHERE;
      this.generateClassicSphere();
    }

    return
  }

  generateClassicSphere() {
    this.positions = new Array(3 * (this.DIV + 1) * (this.DIV + 1));
    this.indices = this.mode == DRAW_MODE.TRIANGLE_STRIP ? new Array(4 * this.DIV * this.DIV) : new Array(3 * this.DIV * this.DIV);
    let iterator = 0;

    // Generate coordinates
    for (let j = 0; j <= this.DIV; ++j) {
      let phi = j * Math.PI / this.DIV;
      let sj = Math.sin(phi);
      let cj = Math.cos(phi);

      for (let i = 0; i <= this.DIV; ++i) {
        let theta = 2 * i * Math.PI / this.DIV;
        let si = Math.sin(theta);
        let ci = Math.cos(theta);

        this.positions[iterator++] = this.center.x + (this.R * ci * sj); // X
        this.positions[iterator++] = this.center.y + (this.R * si * sj); // Y
        this.positions[iterator++] = this.center.z + (this.R * cj); // Z
      }
    }

    iterator = 0;

    // Generate indices (for *_STRIP methods)
    for (let j = 0; j < this.DIV; ++j) {
      for (let i = 0; i < this.DIV; ++i) {
        let p1 = j * (this.DIV + 1) + i;
        let p2 = p1 + (this.DIV + 1);

        this.indices[iterator++] = p1;
        this.indices[iterator++] = p1 + 1;

        if (this.mode == DRAW_MODE.TRIANGLE_STRIP) {
          this.indices[iterator++] = p2;
        }

        this.indices[iterator++] = p2 + 1;
      }
    }
    return;
  }

  generateSpiralSphere() {
    const N = this.DIV;

    this.positions = new Array(3 * N);
    let iterator = 0;

    // Generate coordinates
    let theta = 0;

    for (let k = 1; k <= N; ++k) {
      let h = -1 + (2 * k - 2) / (N - 1);
      let phi = Math.acos(h);

      let s_phi = Math.sin(phi);
      let c_phi = Math.cos(phi);

      theta = (theta + 3.8 / Math.sqrt(N * (1 - h * h))) % (2 * Math.PI);
      if (k == N || k == 1)
        theta = 0;

      let s_theta = Math.sin(theta);
      let c_theta = Math.cos(theta);

      this.positions[iterator++] = this.center.x + (this.R * s_phi * c_theta); // X
      this.positions[iterator++] = this.center.y + (this.R * s_phi * s_theta); // Y
      this.positions[iterator++] = this.center.z + (this.R * c_phi); // Z
    }

    // Generate indices (for POINTS method)
    this.indices = new Array(N);
    for (let i = 0; i < N; i++)
      this.indices[i] = i;

    return;
  }
}


class Camera {
  // TODO handler on keys to change camera pos (wasd, qe)
  constructor() {
    this.pos = DEFAULT_CAMERA_POS;
    this.lookAt = DEFAULT_CAMERA_LOOK_AT;
    this.up = DEFAULT_CAMERA_UP;
  }
}


class Graphic {

  constructor(sphere) { // TODO sphere to interface
    this.canvas = document.getElementById('webgl');

    this.points = null;
    this.isolines = [];
    this.sphere = sphere;
    this.camera = new Camera();

    // params needed for scene rotation by mouse movement
    this.drag = false;
    this.old_x;
    this.old_y;
    this.dX = 0;
    this.dY = 0;
    this.THETA = 0;
    this.PHI = 0;

    this.canvas.onmousedown = (e) => {
      this.drag = true;
      this.old_x = e.pageX;
      this.old_y = e.pageY;
      e.preventDefault();
      return false;
    };
    this.canvas.onmouseup = (e) => {
      this.drag = false;
    };
    this.canvas.onmousemove = (e) => {
      if (!this.drag) return false;
      this.dX = (e.pageX - this.old_x) * 2 * Math.PI / 10;
      this.dY = (e.pageY - this.old_y) * 2 * Math.PI / 10;

      this.THETA += this.dX;
      this.PHI += this.dY;

      this.THETA %= 360;
      this.PHI %= 360;

      this.old_x = e.pageX;
      this.old_y = e.pageY;
      e.preventDefault();
    };
    //

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

    //this.gl.clearDepth(0.5);
    this.gl.enable(this.gl.DEPTH_TEST);
    //this.gl.depthFunc(this.gl.LEQUAL);

    return;
  }

  setupPoints(positions) {
    this.points = {
      positions: positions,
      indices: [...Array(positions.length / 3).keys()],
      color: POINTS_COLOR,
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
    case DRAW_MODE.TRIANGLE_STRIP:
      this.gl.drawElements(this.gl.TRIANGLE_STRIP, vertex_n, this.gl.UNSIGNED_INT, 0);
      break;
    }

    return;
  }

  updateState() {
    //this.camera.updateState();
  }

  tick() {
    this.updateState();
    this.clear();

    this.draw(this.sphere);
    if (this.points) {
      this.draw(this.points);
    }
    if (this.isolines.length) {
      this.isolines.forEach(
        (line) => {this.draw(line)}
      )
    }

    requestAnimationFrame(() => {
      this.tick()
    });

    return;
  }

}


//let sphere = new Sphere(SPIRAL_SPHERE, DRAW_MODE.POINTS, 100000, 1, {x:0, y:0, z:0});
let sphere = new Sphere(CLASSIC_SPHERE, DRAW_MODE.LINE_STRIP, 46, 1, {x:0, y:0, z:0});
let graphic = new Graphic(sphere);

function main() {
  graphic.tick();
}
