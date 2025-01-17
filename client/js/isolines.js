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

    graphic.setupMeanPoint(mean);

    return [mean, covMatrix];
}

function prepareIsolineAPIQuery() {
    let [mean, covMatrix] = getDistributionParams();

    let meanStr = mean.join(',');
    let covStr = covMatrix.join(',');
    let ratios = ISOLINE_AREA_RATIO.join(',')

    return SERVER_URL + ISOLINE_API + `?mean=${meanStr}&cov=${covStr}&ratio=${ratios}`;
}

function getIsolines() {
    let url = prepareIsolineAPIQuery()
    let xhr = new XMLHttpRequest();

    xhr.open("GET", url);
    xhr.onreadystatechange = () => {
        if (xhr.readyState != 4) return;
        if (xhr.status != 200) {
            alert(xhr.status + ': ' + xhr.statusText);
            return;
        }

        let jsonResponse = JSON.parse(xhr.responseText);
        if (jsonResponse.code != 200) {
            alert(jsonResponse.code + ': ' + jsonResponse.body);
        }

        let isolines = []
        jsonResponse.body.forEach(e => isolines.push(e[1]))
        graphic.setupIsolines(isolines);
    }

    xhr.send();
}

function clearIsolines() {
    graphic.setupIsolines();
    graphic.setupDebugIsolines();
}
