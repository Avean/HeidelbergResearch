import QtQuick

// Overlay of the final profiles of all realizations: converged ones in blue,
// not converged ones in red, all drawn with thin translucent lines.
Canvas {
    id: chart

    property var xValues: []
    property var profiles: []
    property var converged: []
    property string title: ""
    property string subtitle: ""

    onProfilesChanged: requestPaint()
    onXValuesChanged: requestPaint()
    onWidthChanged: requestPaint()
    onHeightChanged: requestPaint()

    onPaint: {
        const context = getContext("2d")
        context.reset()
        context.fillStyle = "#ffffff"
        context.fillRect(0, 0, width, height)

        const left = 42
        const right = 12
        const top = subtitle.length > 0 ? 32 : 24
        const bottom = 36
        const plotWidth = Math.max(1, width - left - right)
        const plotHeight = Math.max(1, height - top - bottom)
        const xs = xValues || []
        const lines = profiles || []

        context.fillStyle = "#20252d"
        context.font = "bold 13px sans-serif"
        context.fillText(title, left, 17)
        context.font = "11px sans-serif"
        context.fillStyle = "#68717d"
        if (subtitle.length > 0)
            context.fillText(subtitle, left, 29)

        context.strokeStyle = "#9ba6b5"
        context.lineWidth = 1
        context.beginPath()
        context.moveTo(left, top)
        context.lineTo(left, top + plotHeight)
        context.lineTo(left + plotWidth, top + plotHeight)
        context.stroke()

        if (xs.length < 2 || lines.length === 0)
            return

        let ymin = Infinity
        let ymax = -Infinity
        for (let line = 0; line < lines.length; ++line) {
            const profile = lines[line]
            for (let index = 0; index < profile.length; ++index) {
                const value = profile[index]
                if (value === null || !isFinite(value))
                    continue
                ymin = Math.min(ymin, value)
                ymax = Math.max(ymax, value)
            }
        }
        if (!isFinite(ymin) || !isFinite(ymax))
            return
        if (ymax - ymin < 1e-12) {
            ymin -= 0.5
            ymax += 0.5
        }
        const pad = 0.05 * (ymax - ymin)
        ymin -= pad
        ymax += pad

        const xmin = xs[0]
        const xmax = xs[xs.length - 1]
        const xScale = plotWidth / Math.max(1e-12, xmax - xmin)
        const yScale = plotHeight / (ymax - ymin)

        context.fillStyle = "#68717d"
        context.font = "10px sans-serif"
        context.fillText(Number(ymax).toPrecision(3), 2, top + 4)
        context.fillText(Number(ymin).toPrecision(3), 2, top + plotHeight + 3)

        context.lineWidth = 1
        for (let line = 0; line < lines.length; ++line) {
            const profile = lines[line]
            const ok = line < converged.length ? converged[line] : true
            context.strokeStyle = ok ? "rgba(37, 99, 235, 0.28)" : "rgba(220, 38, 38, 0.55)"
            context.beginPath()
            let drawing = false
            for (let index = 0; index < profile.length && index < xs.length; ++index) {
                const value = profile[index]
                if (value === null || !isFinite(value)) {
                    drawing = false
                    continue
                }
                const x = left + (xs[index] - xmin) * xScale
                const y = top + plotHeight - (value - ymin) * yScale
                if (drawing)
                    context.lineTo(x, y)
                else
                    context.moveTo(x, y)
                drawing = true
            }
            context.stroke()
        }

        context.fillStyle = "#4b5563"
        context.font = "9px sans-serif"
        context.textAlign = "center"
        context.fillText(Number(xmin).toPrecision(4), left, top + plotHeight + 13)
        context.fillText(Number((xmin + xmax) / 2).toPrecision(4), left + plotWidth / 2, top + plotHeight + 13)
        context.fillText(Number(xmax).toPrecision(4), left + plotWidth, top + plotHeight + 13)
        context.textAlign = "start"
    }
}
