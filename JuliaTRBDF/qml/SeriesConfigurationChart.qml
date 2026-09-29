import QtQuick

// Horizontal bars, one per head configuration (e.g. "H------", "---H---"),
// coloured by the number of heads, followed by the not-converged count.
Canvas {
    id: chart

    property var labels: []
    property var counts: []
    property var heads: []
    property int notConverged: 0
    property string title: ""

    readonly property int headerHeight: 28
    // Rows shrink when there are more configurations than fit the height.
    readonly property real rowHeight: Math.max(9, Math.min(20, (height - headerHeight - 10) / Math.max(1, rowCount)))
    readonly property int rowCount: (labels ? labels.length : 0) + 1
    readonly property var barPalette: ["#64748b", "#2563eb", "#7c3aed", "#db2777", "#ea580c", "#16a34a", "#0891b2"]


    onLabelsChanged: requestPaint()
    onCountsChanged: requestPaint()
    onNotConvergedChanged: requestPaint()
    onWidthChanged: requestPaint()
    onHeightChanged: requestPaint()

    onPaint: {
        const context = getContext("2d")
        context.reset()
        context.fillStyle = "#ffffff"
        context.fillRect(0, 0, width, height)

        context.fillStyle = "#20252d"
        context.font = "bold 13px sans-serif"
        context.fillText(title, 42, 17)

        const names = (labels || []).concat(["Not converged"])
        const values = (counts || []).concat([notConverged])
        let longest = 0
        for (let index = 0; index < names.length; ++index)
            longest = Math.max(longest, String(names[index]).length)

        const labelWidth = Math.min(width * 0.45, 12 + longest * 7.3)
        const countWidth = 44
        const barLeft = labelWidth + 8
        const barWidth = Math.max(1, width - barLeft - countWidth)
        let maximum = 1
        for (let index = 0; index < values.length; ++index)
            maximum = Math.max(maximum, Number(values[index]) || 0)

        for (let index = 0; index < names.length; ++index) {
            const y = headerHeight + index * rowHeight
            const value = Number(values[index]) || 0
            const notConvergedRow = index === names.length - 1
            const headCount = notConvergedRow ? -1 : Number(heads[index])

            context.fillStyle = notConvergedRow ? "#9f252b" : "#20252d"
            context.font = notConvergedRow
                            ? Math.round(Math.min(11, rowHeight - 4)) + "px sans-serif"
                            : Math.round(Math.min(12, rowHeight - 3)) + "px monospace"
            context.textBaseline = "middle"
            context.fillText(String(names[index]), 6, y + rowHeight / 2)

            context.fillStyle = notConvergedRow
                                ? "#dc2626"
                                : barPalette[Math.min(headCount, barPalette.length - 1)]
            context.fillRect(barLeft, y + 2, barWidth * value / maximum, Math.max(2, rowHeight - 5))

            context.fillStyle = "#4b5563"
            context.font = "11px sans-serif"
            context.fillText(String(value), barLeft + barWidth * value / maximum + 6, y + rowHeight / 2)
        }
        context.textBaseline = "alphabetic"
    }
}
