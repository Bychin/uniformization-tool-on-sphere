function prepareStatsAPIJsonData() {
    let [mean, covMatrix] = getDistributionParams();

    let meanStr = mean.join(',');
    let covStr = covMatrix.join(',');
    let points = graphic.points.positions.join(',') // TODO empty checks

    return JSON.stringify({
        points: points,
        mean: meanStr,
        cov: covStr,
    });
}


function getStats() {
    let xhr = new XMLHttpRequest();
    let data = prepareStatsAPIJsonData()

    xhr.open("POST", SERVER_URL + STATS_API)
    xhr.setRequestHeader('Content-type', 'application/json');
    xhr.onreadystatechange = () => {
        if (xhr.readyState != 4) return;
        if (xhr.status != 200) {
            alert(xhr.status + ': ' + xhr.statusText);
            return;
        }

        let jsonResponse = JSON.parse(xhr.responseText);
        console.log(jsonResponse)

        let isolines = []
        for (let sp of jsonResponse.s) {
            isolines.push(sp.isoline.flat());
        }
        graphic.setupIsolines(isolines);

        let points = []
        for (let sp of jsonResponse.s) {
            points.push(sp.points.flat())
        };
        graphic.setupMorePoints(points.flat());
    };

    xhr.send(data);
}
