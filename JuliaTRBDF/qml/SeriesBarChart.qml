import QtQuick

Canvas {
    id: chart

    property var values: []
    property var labels: []
    property string title: ""
    property string subtitle: ""
    property real xMinimum: NaN
    property real xMaximum: NaN
    property color barColor: "#3b82f6"

    onValuesChanged: requestPaint()
    onLabelsChanged: requestPaint()
    onWidthChanged: requestPaint()
    onHeightChanged: requestPaint()

    onPaint: {
        const context = getContext("2d")
        context.reset()
        context.fillStyle = "#ffffff"
        context.fillRect(0, 0, width, height)

        const left = 42
        const right = 12
        const top = 32
        const bottom = labels.length > 0 ? 36 : (isFinite(xMinimum) && isFinite(xMaximum) ? 36 : 24)
        const plotWidth = Math.max(1, width - left - right)
        const plotHeight = Math.max(1, height - top - bottom)
        const data = values || []
        let maximum = 0
        for (let index = 0; index < data.length; ++index)
            maximum = Math.max(maximum, Number(data[index]) || 0)
        maximum = Math.max(1, maximum)

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

        context.fillStyle = "#68717d"
        context.font = "10px sans-serif"
        context.fillText(String(maximum), 4, top + 4)
        context.fillText("0", 23, top + plotHeight + 3)

        if (data.length === 0)
            return

        const gap = data.length > 32 ? 0 : 2
        const barWidth = Math.max(1, (plotWidth - gap * (data.length - 1)) / data.length)
        context.fillStyle = barColor

        for (let index = 0; index < data.length; ++index) {
            const value = Number(data[index]) || 0
            const barHeight = plotHeight * value / maximum
            const x = left + index * (barWidth + gap)
            context.fillRect(x, top + plotHeight - barHeight, barWidth, barHeight)
        }

        if (labels.length > 0 && data.length <= 14) {
            context.fillStyle = "#4b5563"
            context.font = "9px sans-serif"
            context.textAlign = "center"
            for (let index = 0; index < labels.length; ++index) {
                const x = left + (index + 0.5) * (plotWidth / labels.length)
                context.save()
                context.translate(x, top + plotHeight + 11)
                context.rotate(-Math.PI / 5)
                context.fillText(String(labels[index]), 0, 0)
                context.restore()
            }
            context.textAlign = "start"
        } else if (isFinite(xMinimum) && isFinite(xMaximum)) {
            context.fillStyle = "#4b5563"
            context.font = "9px sans-serif"
            context.textAlign = "center"
            context.fillText(Number(xMinimum).toPrecision(4), left, top + plotHeight + 13)
            context.fillText(Number((xMinimum + xMaximum) / 2).toPrecision(4), left + plotWidth / 2, top + plotHeight + 13)
            context.fillText(Number(xMaximum).toPrecision(4), left + plotWidth, top + plotHeight + 13)
            context.textAlign = "start"
        }
    }
}
