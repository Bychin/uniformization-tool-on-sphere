function uploadPointsFile() {
    if (!window.FileReader) {
        alert("Your browser is not supported")
    }

    let input = document.getElementById("pointsFile");
    let reader = new FileReader();
    if (input.files.length) {
        let textFile = input.files[0];
        
        reader.onload = (e) => {processPointsFile(e)};
        reader.onerror = (e) => {alert("Error while reading points file")};

        reader.readAsText(textFile);
    } else {
        alert("Upload a file before continuing");
    }
};

function processPointsFile(e) {
    let file = e.target.result;
    if (!file || !file.length) {
        return;
    }

    let positions = [];
    let pointsRaw = file.split("\n");
    for (let point of pointsRaw) {
        let coords = point.split(" "); // get X, Y, Z coordinates of a point
        let normalizedCoords = normalizePoint(coords);

        for (let coord of normalizedCoords) {
            positions.push(coord);
        }
    }

    graphic.setupPoints(positions);
}

// normalizePoint normalizes point in 3 dimensional space with X, Y, Z coordinates
// which must be given as Array(3)
function normalizePoint(coords) {
    let length = Math.sqrt(coords[0] * coords[0] + coords[1] * coords[1] + coords[2] * coords[2]);

    for (let i = 0; i < 3; ++i) {
        coords[i] /= length;
    }

    return coords;
}
