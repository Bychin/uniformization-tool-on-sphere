const SERVER_URL = "http://127.0.0.1:8080"
const ISOLINE_API = "/api/isoline"

const ISOLINE_AREA_RATIO = [0.9, 0.7, 0.5, 0.3, 0.1];


// "1/5" -> "0.2"
function convertFraction(numRaw) {
    if (!numRaw.includes('/')) {
        return numRaw;
    }

    let [n, d] = numRaw.split('/');
    return parseFloat(n / d);
}


function getDistributionParams() {
    let mean = Array(3); // 3-dimensional vector
    mean[0] = document.getElementById("meanX").value;
    mean[1] = document.getElementById("meanY").value;
    mean[2] = document.getElementById("meanZ").value;
    mean.forEach((c, i, arr) => {arr[i] = convertFraction(c)});

    // symmetrical covariance matrix (stored by lines as upper triangular matrix)
    let covMatrix = Array(6);
    covMatrix[0] = document.getElementById("cov11").value;
    covMatrix[1] = document.getElementById("cov12").value;
    covMatrix[2] = document.getElementById("cov13").value;
    covMatrix[3] = document.getElementById("cov22").value;
    covMatrix[4] = document.getElementById("cov23").value;
    covMatrix[5] = document.getElementById("cov33").value;
    covMatrix.forEach((c, i, arr) => {arr[i] = convertFraction(c)});

    return [mean, covMatrix];
}

function prepareIsolineAPIQuery(ratio) {
    let [mean, covMatrix] = getDistributionParams();

    let meanStr = mean.join(',');
    let covStr = covMatrix.join(',');
    let ratios = ISOLINE_AREA_RATIO.join(',')

    return SERVER_URL + ISOLINE_API + `?mean=${meanStr}&cov=${covStr}&ratio=${ratios}`;
}

function uploadAGD() {
    let url = prepareIsolineAPIQuery(ISOLINE_AREA_RATIO[0])
    let xhr = new XMLHttpRequest();

    xhr.open("GET", url);
    xhr.send();
    xhr.onreadystatechange = () => {
        if (xhr.readyState != 4) return;
        if (xhr.status != 200) {
            alert(xhr.status + ': ' + xhr.statusText);
            return;
        }

        let jsonResponse = JSON.parse(xhr.responseText);
        console.log(jsonResponse)
        let isolinesDict = jsonResponse.isolines;
        let isolines = []
        for (let ratio in isolinesDict) {
            isolines.push(isolinesDict[ratio].flat());
        }
        graphic.setupIsolines(isolines);
    }
}
