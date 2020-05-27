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
        let resp = jsonResponse.body;

        printTestsResults(resp);

        if (resp.s_stats_debug_info === undefined)
            return;

        drawDebugElements(resp);
    };

    xhr.send(data);
}

function printTestsResults(response) {
    let area = document.getElementById("results-area");
    area.value = `Nonparametric tests results:
t-statistic:
KS-test: Measure: ${(response.t_results.ks_measure || NaN).toFixed(8)}, Estimate: ${(response.t_results.ks_est || NaN).toFixed(8)}
AD-test: Measure: ${(response.t_results.ad_measure || NaN).toFixed(8)}, Estimate: ${(response.t_results.ad_est || NaN).toFixed(8)}
s-statistic:
KS-test: Measure: ${(response.s_results.ks_measure || NaN).toFixed(8)}, Estimate: ${(response.s_results.ks_est || NaN).toFixed(8)}
AD-test: Measure: ${(response.s_results.ad_measure || NaN).toFixed(8)}, Estimate: ${(response.s_results.ad_est || NaN).toFixed(8)}`;
}

function drawDebugElements(response) {
    let debugIsolines = [];
    let debugPoints = [];
    let debugDirPoints = [];

    for (let s_debug of response.s_stats_debug_info) {
        debugIsolines.push(s_debug.isoline);
        debugPoints.push(s_debug.isoline[s_debug.intersection_point].flat());

        let dirPoint = s_debug.clockwise_direction ? s_debug.isoline[s_debug.intersection_point+3] : s_debug.isoline[s_debug.intersection_point-3]; // TODO to cfg?
        if (dirPoint === undefined) {
            if (s_debug.clockwise_direction)
                dirPoint = s_debug.isoline[0];
            else
                dirPoint = s_debug.isoline[s_debug.isoline.length-1];
        }

        debugDirPoints.push(dirPoint.flat());
    }

    graphic.setupDebugIsolines(debugIsolines);
    graphic.setupDebugIntPoints(debugPoints.flat());
    graphic.setupDebugDirPoints(debugDirPoints.flat());
}
