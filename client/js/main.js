let graphic = new Graphic(sphere);

function main() {
    let sphere = new Sphere(CLASSIC_SPHERE, DRAW_MODE.LINE_STRIP, 36, 1, {x:0, y:0, z:0});

    graphic.setupCoordLines();
    graphic.setupSphere(sphere);

    graphic.tick();
}
