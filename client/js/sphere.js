class Sphere {
  constructor(type, mode, div, radius, center) {
    this.mode = mode;
    this.color = COLORS.LIGHT_RED;

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
  }
}
