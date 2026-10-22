/* Single-page behaviour: submit the form in place, show results below, draw charts. */
(function () {
    const form = document.getElementById("analyze-form");
    const box = document.getElementById("query");
    const go = document.getElementById("go");
    const results = document.getElementById("results");
    const page = document.querySelector(".page");

    // If the page came from a plain form post (/analyze), point the address bar
    // back at "/" so a reload or a restored tab starts empty.
    if (window.location.pathname !== "/") history.replaceState(null, "", "/");

    // ---------- Settings pop-out ----------
    const settings = document.querySelector(".settings");
    document.addEventListener("click", event => {
        if (settings.open && !settings.contains(event.target)) settings.open = false;
    });
    document.addEventListener("keydown", event => {
        if (event.key === "Escape") settings.open = false;
    });
    const reset = settings.querySelector("[data-reset]");
    reset.addEventListener("click", () => {
        ["window", "threshold", "pattern", "allowed"].forEach(name => {
            form.elements[name].value = reset.dataset[name];
        });
    });

    // ---------- Input ----------
    const grow = () => { box.style.height = "auto"; box.style.height = Math.min(box.scrollHeight, 132) + "px"; };
    box.addEventListener("input", grow);
    grow();

    const chips = document.querySelectorAll(".chip");
    let activeChip = null;
    chips.forEach(chip => chip.addEventListener("click", () => {
        if (busy) return;
        box.value = chip.hasAttribute("data-fill-example") ? window.EXAMPLE_SEQUENCE : chip.dataset.fill;
        grow();
        activeChip = chip;
        form.requestSubmit();
    }));
    // Enter submits; Shift+Enter adds a new line.
    box.addEventListener("keydown", event => {
        if (event.key === "Enter" && !event.shiftKey) { event.preventDefault(); form.requestSubmit(); }
    });

    // Same rule as the server: one short token containing a digit is an accession.
    const looksLikeAccession = text => /^[A-Za-z0-9_.]{2,30}$/.test(text) && /\d/.test(text);

    function showLoader() {
        const query = box.value.trim();
        const residues = query.replace(/^>.*$/gm, "").replace(/[^A-Za-z]/g, "").length;
        const title = looksLikeAccession(query)
            ? `Fetching <span class="mono">${query.toUpperCase().replace(/[^A-Z0-9_.]/g, "")}</span> from NCBI`
            : `Analyzing a ${residues}-residue sequence`;
        results.innerHTML = `
            <div class="loader" role="status">
                <div class="loader-head">
                    <span class="loader-spin" aria-hidden="true"></span>
                    <div>
                        <p class="loader-title">${title}…</p>
                        <p class="loader-sub">Scoring hydrophobicity and searching for the motif</p>
                    </div>
                </div>
                <div class="skeleton">
                    <span style="width: 42%"></span><span style="width: 100%; height: 120px"></span>
                    <span style="width: 70%"></span><span style="width: 86%"></span>
                </div>
            </div>`;
    }

    let busy = false;
    form.addEventListener("submit", async event => {
        event.preventDefault();
        if (busy) return;
        busy = true;
        page.classList.remove("is-idle");   // the search glides from the centre to the top
        go.classList.add("is-loading");
        go.disabled = true;
        chips.forEach(chip => chip.classList.toggle("is-active", chip === activeChip));
        showLoader();
        try {
            const response = await fetch(form.action, {
                method: "POST",
                body: new FormData(form),
                headers: { "X-Requested-With": "fetch" },
            });
            results.innerHTML = await response.text();
        } catch (error) {
            results.innerHTML = '<p class="error" role="alert">Could not reach the server. Check your connection and try again.</p>';
        } finally {
            busy = false;
            activeChip = null;
            go.classList.remove("is-loading");
            go.disabled = false;
            chips.forEach(chip => chip.classList.remove("is-active"));
        }
        setupResults();
        window.scrollTo({ top: 0, behavior: "smooth" });
    });

    // ---------- Charts ----------
    const narrow = () => window.matchMedia("(max-width: 640px)").matches;
    const css = getComputedStyle(document.documentElement);
    const color = name => css.getPropertyValue(name).trim();
    const C = {
        ink: color("--ink"), muted: color("--muted"), grid: color("--grid"), line: color("--line"),
        red: color("--red"), redFill: color("--red-fill"),
        blue: color("--blue"), blueBand: color("--blue-band"), surface: color("--surface"),
    };
    if (window.Chart) {
        Chart.defaults.font.family = getComputedStyle(document.body).fontFamily;
        Chart.defaults.font.size = 12;
        Chart.defaults.color = C.muted;
    }

    // White background so downloaded PNGs are not transparent.
    const background = {
        id: "background",
        beforeDraw(chart) {
            chart.ctx.save();
            chart.ctx.fillStyle = C.surface;
            chart.ctx.fillRect(0, 0, chart.width, chart.height);
            chart.ctx.restore();
        },
    };
    // Shaded bands marking qualifying motifs on the profile.
    const motifBands = {
        id: "motifBands",
        beforeDatasetsDraw(chart, _args, options) {
            const { ctx, chartArea, scales } = chart;
            ctx.save();
            let lastLabelEnd = -Infinity;
            for (const band of options.bands || []) {
                const x1 = scales.x.getPixelForValue(band.start - 0.5);
                const width = Math.max(scales.x.getPixelForValue(band.end + 0.5) - x1, 3);
                ctx.fillStyle = C.blueBand;
                ctx.fillRect(x1, chartArea.top, width, chartArea.bottom - chartArea.top);
                // Skip a name that would overlap the previous one; the band itself still shows.
                ctx.font = '600 11px "JetBrains Mono", monospace';
                const labelWidth = ctx.measureText(band.motif).width;
                if (x1 + width / 2 - labelWidth / 2 < lastLabelEnd + 6) continue;
                lastLabelEnd = x1 + width / 2 + labelWidth / 2;
                ctx.fillStyle = C.blue;
                ctx.textAlign = x1 > chartArea.right - 60 ? "right" : "center";
                const x = ctx.textAlign === "right" ? x1 + width : x1 + width / 2;
                ctx.fillText(band.motif, x, chartArea.top - 6);
            }
            ctx.restore();
        },
    };
    // Shades the motif's own residues on a per-motif plot (category axis).
    const indexBand = {
        id: "indexBand",
        beforeDatasetsDraw(chart, _args, options) {
            if (options.from == null) return;
            const { ctx, chartArea, scales } = chart;
            const step = scales.x.getPixelForValue(1) - scales.x.getPixelForValue(0);
            const x1 = scales.x.getPixelForValue(options.from) - step / 2;
            const x2 = scales.x.getPixelForValue(options.to) + step / 2;
            ctx.save();
            ctx.fillStyle = C.blueBand;
            ctx.fillRect(Math.max(x1, chartArea.left), chartArea.top, Math.min(x2, chartArea.right) - Math.max(x1, chartArea.left), chartArea.bottom - chartArea.top);
            ctx.restore();
        },
    };


    function profileChart(data, sequence) {
        const canvas = document.getElementById("profile-chart");
        if (!canvas || !data.profile.length) return;
        const half = Math.floor((data.window - 1) / 2);
        new Chart(canvas, {
            type: "line",
            data: {
                datasets: [
                    {
                        data: data.profile, borderColor: C.ink, borderWidth: 1.4,
                        pointRadius: 0, pointHoverRadius: 4, tension: 0,
                        fill: { target: { value: data.threshold }, above: C.redFill, below: "transparent" },
                    },
                    {
                        data: [{ x: 1, y: data.threshold }, { x: data.length, y: data.threshold }],
                        borderColor: C.red, borderDash: [5, 4], borderWidth: 1.2, pointRadius: 0, pointHoverRadius: 0,
                    },
                ],
            },
            options: {
                maintainAspectRatio: false, animation: false, parsing: false,
                layout: { padding: { top: 16 } },
                interaction: { mode: "nearest", axis: "x", intersect: false },
                scales: {
                    x: { type: "linear", min: 1, max: data.length, grid: { display: false }, border: { color: C.line },
                         title: { display: true, text: "Residue position" }, ticks: { maxTicksLimit: narrow() ? 6 : 12 } },
                    y: { suggestedMin: -3, suggestedMax: 3, grid: { color: C.grid }, border: { display: false },
                         title: { display: !narrow(), text: "Mean KD score" } },
                },
                plugins: {
                    legend: { display: false },
                    tooltip: {
                        filter: item => item.datasetIndex === 0,
                        callbacks: {
                            title: items => {
                                const pos = items[0].parsed.x;
                                return `${sequence[pos - 1]}${pos}, window ${pos - half}–${pos - half + data.window - 1}`;
                            },
                            label: item => `Score ${item.parsed.y.toFixed(2)}`,
                        },
                    },
                    motifBands: { bands: data.motifs.filter(m => m.qualifying) },
                },
            },
            plugins: [background, motifBands],
        });
    }

    function motifCharts(data) {
        data.motifs.forEach((motif, index) => {
            const canvas = document.getElementById(`motif-chart-${index}`);
            if (!canvas) return;
            const pass = motif.scores.map(score => score >= data.threshold - 1e-9);
            const inMotif = motif.labels.map(label => {
                const pos = parseInt(label.match(/\((\d+)\)/)[1], 10);
                return pos >= motif.start && pos <= motif.end;
            });
            const first = inMotif.indexOf(true);
            const last = inMotif.lastIndexOf(true);
            const lo = Math.min(...motif.scores, data.threshold);
            const hi = Math.max(...motif.scores, data.threshold);
            const pad = Math.max((hi - lo) * 0.15, 0.2);
            new Chart(canvas, {
                type: "line",
                data: {
                    labels: motif.labels,
                    datasets: [
                        {
                            data: motif.scores, borderColor: C.ink, borderWidth: 1.5, tension: 0,
                            pointRadius: 4, pointHoverRadius: 6, pointBorderWidth: 1.5,
                            pointBackgroundColor: pass.map(p => (p ? C.red : C.surface)),
                            pointBorderColor: pass.map(p => (p ? C.red : C.ink)),
                        },
                        {
                            data: motif.labels.map(() => data.threshold),
                            borderColor: C.red, borderDash: [5, 4], borderWidth: 1.2, pointRadius: 0, pointHoverRadius: 0,
                        },
                    ],
                },
                options: {
                    maintainAspectRatio: false, animation: false,
                    interaction: { mode: "index", intersect: false },
                    scales: {
                        x: {
                            grid: { display: false }, border: { color: C.line },
                            ticks: { autoSkip: false, maxRotation: 60, font: { size: narrow() ? 9.5 : 11 },
                                     color: ctx => (inMotif[ctx.index] ? C.blue : C.muted) },
                        },
                        y: {
                            min: Math.floor((lo - pad) * 10) / 10, max: Math.ceil((hi + pad) * 10) / 10,
                            grid: { color: C.grid }, border: { display: false }, title: { display: !narrow(), text: "Mean KD score" },
                        },
                    },
                    plugins: {
                        legend: { display: false },
                        tooltip: {
                            filter: item => item.datasetIndex === 0,
                            callbacks: {
                                title: items => `Window ${motif.ranges[items[0].dataIndex]}`,
                                label: item => `Score ${item.parsed.y.toFixed(2)}`,
                            },
                        },
                        indexBand: first === -1 ? {} : { from: first, to: last },
                    },
                },
                plugins: [background, indexBand],
            });
        });
    }

    // ---------- Result actions ----------
    function copyText(text, button) {
        const original = button.textContent;
        const done = () => { button.textContent = "Copied"; setTimeout(() => (button.textContent = original), 1500); };
        if (navigator.clipboard && window.isSecureContext) {
            navigator.clipboard.writeText(text).then(done);
        } else {
            const area = document.createElement("textarea");
            area.value = text;
            document.body.appendChild(area);
            area.select();
            document.execCommand("copy");
            area.remove();
            done();
        }
    }

    function setupResults() {
        if (window.Chart) Object.values(Chart.instances).forEach(chart => chart.destroy());
        const node = results.querySelector(".analysis-data");
        if (!node) return;
        const payload = JSON.parse(node.textContent);
        if (window.Chart) {
            profileChart(payload.chart, payload.sequence);
            motifCharts(payload.chart);
        }

        const fasta = results.querySelector("[data-copy-fasta]");
        if (fasta) fasta.addEventListener("click", () =>
            copyText(`>${payload.fasta_header}\n${payload.sequence.match(/.{1,60}/g).join("\n")}\n`, fasta));
        results.querySelectorAll("[data-download]").forEach(button => button.addEventListener("click", () => {
            const chart = window.Chart && Chart.getChart(button.dataset.download);
            if (!chart) return;
            const a = document.createElement("a");
            a.href = chart.toBase64Image("image/png", 1);
            a.download = `${button.dataset.name}.png`;
            a.click();
        }));
    }

    setupResults();
})();
