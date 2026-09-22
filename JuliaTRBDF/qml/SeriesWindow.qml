import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import QtQuick.Window
import jlqml

// Series mode window: series setup on the left, head statistics on the right.
// It is visible exactly while ui.seriesMode is on; closing it leaves series
// mode, which is refused while a series is running.
Window {
    id: seriesWindow

    property var variables: JSON.parse(ui.variablesJson)
    property var perturbations: JSON.parse(ui.seriesPerturbationsJson)
    property var results: {
        try {
            return JSON.parse(ui.seriesResultsJson)
        } catch (error) {
            return []
        }
    }
    property var panelNames: Array.from({length: ui.segmentCount}, (_, index) => String(index + 1))

    visible: ui.seriesMode
    width: 1400
    height: 820
    minimumWidth: 900
    minimumHeight: 600
    title: ui.seriesRunning
           ? "Series — " + ui.seriesCompletedRuns + " / " + ui.seriesTotalRuns
           : "Series"
    color: "#eef1f5"

    onClosing: function(close) {
        close.accepted = false
        if (!ui.seriesRunning)
            Julia.setSeriesMode(false)
    }

    ColumnLayout {
        anchors.fill: parent
        spacing: 0

        Rectangle {
            Layout.fillWidth: true
            Layout.preferredHeight: 50
            color: "#20252d"

            RowLayout {
                anchors.fill: parent
                anchors.leftMargin: 14
                anchors.rightMargin: 8

                Label {
                    text: "Simulation series"
                    color: "white"
                    font.bold: true
                    font.pixelSize: 16
                }

                Label {
                    Layout.fillWidth: true
                    Layout.leftMargin: 12
                    text: ui.seriesStatus
                    color: ui.seriesRunning ? "#93c5fd" : "#cbd3dc"
                    elide: Text.ElideRight
                }

                ToolButton {
                    text: "Close series mode"
                    enabled: !ui.seriesRunning
                    palette.buttonText: "white"
                    onClicked: Julia.setSeriesMode(false)
                    ToolTip.visible: hovered && ui.seriesRunning
                    ToolTip.text: "Stop the series first"
                }
            }
        }

        RowLayout {
            Layout.fillWidth: true
            Layout.fillHeight: true
            spacing: 0

            // -------------------------------------------------- setup
            Rectangle {
                Layout.preferredWidth: 540
                Layout.fillHeight: true
                color: "#f3f5f8"
                border.color: "#c4ccd7"
                border.width: 1

                ScrollView {
                    anchors.fill: parent
                    contentWidth: availableWidth
                    clip: true

                    ColumnLayout {
                        width: parent.width
                        spacing: 9

                        // A plain rectangle rather than a styled Button: the
                        // native Windows style ignores custom backgrounds.
                        Rectangle {
                            id: startStopButton
                            property bool active: ui.seriesRunning || seriesWindow.perturbations.length > 0
                            Layout.leftMargin: 10
                            Layout.rightMargin: 10
                            Layout.topMargin: 12
                            Layout.fillWidth: true
                            Layout.preferredHeight: 56
                            radius: 7
                            color: !active
                                   ? "#9ca3af"
                                   : ui.seriesRunning
                                     ? (startStopMouse.pressed ? "#991b1b" : startStopMouse.containsMouse ? "#b91c1c" : "#dc2626")
                                     : (startStopMouse.pressed ? "#14532d" : startStopMouse.containsMouse ? "#15803d" : "#16a34a")

                            Text {
                                anchors.centerIn: parent
                                text: !ui.seriesRunning
                                      ? "Start series"
                                      : ui.seriesSingleRun
                                        ? "Stop run one"
                                        : "Stop series   (" + ui.seriesCompletedRuns + " / " + ui.seriesTotalRuns + ")"
                                color: "white"
                                font.bold: true
                                font.pixelSize: 18
                            }

                            MouseArea {
                                id: startStopMouse
                                anchors.fill: parent
                                enabled: startStopButton.active
                                hoverEnabled: true
                                cursorShape: Qt.PointingHandCursor
                                onClicked: ui.seriesRunning ? Julia.stopSeries() : Julia.startSeries()
                            }
                        }

                        // Diagnostics: one more realization, added to the
                        // statistics, whose final state stays on the main plots.
                        Rectangle {
                            id: runOneButton
                            property bool active: !ui.seriesRunning && seriesWindow.perturbations.length > 0
                            Layout.leftMargin: 10
                            Layout.rightMargin: 10
                            Layout.fillWidth: true
                            Layout.preferredHeight: 40
                            radius: 7
                            color: !active
                                   ? "#c7ccd4"
                                   : runOneMouse.pressed ? "#3730a3" : runOneMouse.containsMouse ? "#4338ca" : "#4f46e5"

                            Text {
                                anchors.centerIn: parent
                                text: "Run one"
                                color: "white"
                                font.bold: true
                                font.pixelSize: 15
                            }

                            MouseArea {
                                id: runOneMouse
                                anchors.fill: parent
                                enabled: runOneButton.active
                                hoverEnabled: true
                                cursorShape: Qt.PointingHandCursor
                                onClicked: Julia.runOneSeries()
                            }

                            ToolTip.visible: runOneMouse.containsMouse
                            ToolTip.text: "Run one realization, add it to the statistics and keep its final state on the main plots"
                        }

                        Label {
                            Layout.leftMargin: 14
                            Layout.rightMargin: 14
                            Layout.fillWidth: true
                            text: !ui.seriesRunning && seriesWindow.perturbations.length === 0
                                  ? "Add at least one perturbation to start a series."
                                  : ui.seriesStatus
                            color: ui.seriesRunning ? "#2563eb" : "#4b5563"
                            wrapMode: Text.WordWrap
                        }

                        ControlSection {
                            title: "New perturbation"
                            Layout.leftMargin: 10
                            Layout.rightMargin: 10
                            enabled: !ui.seriesRunning

                            RowLayout {
                                Layout.fillWidth: true

                                Label { text: "Panel" }

                                ComboBox {
                                    Layout.preferredWidth: 80
                                    model: seriesWindow.panelNames
                                    currentIndex: Math.max(0, ui.seriesSelectedSegment - 1)
                                    onActivated: Julia.selectSeriesSegment(currentIndex + 1)
                                }

                                Label { text: "Variable" }

                                ComboBox {
                                    Layout.fillWidth: true
                                    model: seriesWindow.variables
                                    currentIndex: Math.max(0, ui.seriesSelectedVariable - 1)
                                    onActivated: Julia.selectSeriesVariable(currentIndex + 1)
                                }
                            }

                            RowLayout {
                                Layout.fillWidth: true

                                Label {
                                    text: "Position: " + Number(ui.seriesPosition).toPrecision(5)
                                    Layout.preferredWidth: 138
                                }

                                Slider {
                                    Layout.fillWidth: true
                                    from: 0
                                    to: Math.max(0.000001, Number(ui.seriesSelectedPanelLength))
                                    value: Number(ui.seriesPosition)
                                    onMoved: Julia.setSeriesPosition(value)
                                }
                            }

                            Button {
                                Layout.alignment: Qt.AlignRight
                                text: "Add perturbation"
                                onClicked: Julia.addSeriesPerturbation()
                            }
                        }

                        ControlSection {
                            title: "Series settings"
                            expanded: false
                            Layout.leftMargin: 10
                            Layout.rightMargin: 10
                            enabled: !ui.seriesRunning

                            GridLayout {
                                Layout.fillWidth: true
                                columns: 2
                                columnSpacing: 8
                                rowSpacing: 7

                                Label { text: "Number of runs" }
                                TextField {
                                    selectByMouse: true
                                    validator: IntValidator { bottom: 1 }
                                    text: String(ui.seriesRunCount)
                                    onEditingFinished: Julia.setSeriesRunCount(text)
                                }

                                Label { text: "Maximum time / panel" }
                                TextField {
                                    selectByMouse: true
                                    validator: DoubleValidator { bottom: 0.0000000001; notation: DoubleValidator.ScientificNotation }
                                    text: Number(ui.seriesMaximumTime).toString()
                                    onEditingFinished: Julia.setSeriesMaximumTime(text)
                                }

                                Label { text: "Check interval (time)" }
                                TextField {
                                    selectByMouse: true
                                    validator: DoubleValidator { bottom: 0.0000000001; notation: DoubleValidator.ScientificNotation }
                                    text: Number(ui.seriesCheckInterval).toString()
                                    onEditingFinished: Julia.setSeriesCheckInterval(text)
                                }

                                Label { text: "Tolerance on |du/dt|" }
                                TextField {
                                    selectByMouse: true
                                    validator: DoubleValidator { bottom: 0.000000000000000001; notation: DoubleValidator.ScientificNotation }
                                    text: Number(ui.seriesTolerance).toExponential()
                                    onEditingFinished: Julia.setSeriesTolerance(text)
                                }

                                Label { text: "Consecutive passed checks" }
                                TextField {
                                    selectByMouse: true
                                    validator: IntValidator { bottom: 1 }
                                    text: String(ui.seriesRequiredChecks)
                                    onEditingFinished: Julia.setSeriesRequiredChecks(text)
                                }

                                Label { text: "Series dtmax" }
                                TextField {
                                    selectByMouse: true
                                    validator: DoubleValidator { bottom: 0.0000000001; notation: DoubleValidator.ScientificNotation }
                                    text: Number(ui.seriesDtmax).toString()
                                    onEditingFinished: Julia.setSeriesDtmax(text)
                                }

                                Label { text: "Safety step limit / panel" }
                                TextField {
                                    selectByMouse: true
                                    validator: IntValidator { bottom: 1 }
                                    text: String(ui.seriesMaximumSteps)
                                    onEditingFinished: Julia.setSeriesMaximumSteps(text)
                                }

                                Label { text: "Random seed" }
                                TextField {
                                    selectByMouse: true
                                    text: ui.seriesSeed
                                    onEditingFinished: Julia.setSeriesSeed(text)
                                }

                                Label { text: "Head detection variable" }
                                ComboBox {
                                    model: seriesWindow.variables
                                    currentIndex: Math.max(0, ui.seriesHeadVariable - 1)
                                    onActivated: Julia.setSeriesHeadVariable(currentIndex + 1)
                                }
                            }

                            Switch {
                                Layout.topMargin: 5
                                text: "Live simulation preview"
                                checked: ui.seriesLivePreview
                                onToggled: Julia.setSeriesLivePreview(checked)
                            }
                        }

                        ControlSection {
                            title: "Perturbations"
                            Layout.leftMargin: 10
                            Layout.rightMargin: 10
                            Layout.bottomMargin: 12

                            Label {
                                Layout.fillWidth: true
                                visible: seriesWindow.perturbations.length === 0
                                text: "No perturbations yet."
                                color: "#68717d"
                            }

                            Repeater {
                                model: seriesWindow.perturbations

                                Rectangle {
                                    required property var modelData
                                    Layout.fillWidth: true
                                    implicitHeight: perturbationCard.implicitHeight + 16
                                    color: "#e6eaf0"
                                    radius: 5
                                    border.color: "#bdc6d2"
                                    border.width: 1

                                    ColumnLayout {
                                        id: perturbationCard
                                        anchors.fill: parent
                                        anchors.margins: 8
                                        spacing: 5

                                        RowLayout {
                                            Layout.fillWidth: true

                                            Label {
                                                Layout.fillWidth: true
                                                text: "Perturbation " + modelData.id + "   ·   x = " + Number(modelData.position).toPrecision(5)
                                                font.bold: true
                                            }

                                            ToolButton {
                                                enabled: !ui.seriesRunning
                                                text: "Select"
                                                onClicked: Julia.selectSeriesPerturbation(modelData.id)
                                            }

                                            ToolButton {
                                                enabled: !ui.seriesRunning
                                                text: "Delete"
                                                onClicked: Julia.deleteSeriesPerturbation(modelData.id)
                                            }
                                        }

                                        RowLayout {
                                            Layout.fillWidth: true

                                            Label { text: "Panel" }
                                            ComboBox {
                                                Layout.preferredWidth: 75
                                                enabled: !ui.seriesRunning
                                                model: seriesWindow.panelNames
                                                currentIndex: Math.max(0, modelData.panel - 1)
                                                onActivated: Julia.updateSeriesPerturbation(modelData.id, "panel", currentIndex + 1)
                                            }

                                            Label { text: "Variable" }
                                            ComboBox {
                                                Layout.fillWidth: true
                                                enabled: !ui.seriesRunning
                                                model: seriesWindow.variables
                                                currentIndex: Math.max(0, modelData.variable - 1)
                                                onActivated: Julia.updateSeriesPerturbation(modelData.id, "variable", currentIndex + 1)
                                            }
                                        }

                                        RowLayout {
                                            Layout.fillWidth: true

                                            Label { text: "Position" }

                                            TextField {
                                                Layout.fillWidth: true
                                                selectByMouse: true
                                                enabled: !ui.seriesRunning
                                                text: Number(modelData.position).toPrecision(5)
                                                onEditingFinished: Julia.updateSeriesPerturbation(modelData.id, "position", text)
                                            }
                                        }

                                        GridLayout {
                                            Layout.fillWidth: true
                                            columns: 4
                                            columnSpacing: 6

                                            Label { text: "Width min" }
                                            TextField {
                                                selectByMouse: true
                                                enabled: !ui.seriesRunning
                                                text: Number(modelData.widthMin).toFixed(2)
                                                onEditingFinished: Julia.updateSeriesPerturbation(modelData.id, "widthMin", text)
                                            }
                                            Label { text: "Width max" }
                                            TextField {
                                                selectByMouse: true
                                                enabled: !ui.seriesRunning
                                                text: Number(modelData.widthMax).toFixed(2)
                                                onEditingFinished: Julia.updateSeriesPerturbation(modelData.id, "widthMax", text)
                                            }

                                            Label { text: "Height min" }
                                            TextField {
                                                selectByMouse: true
                                                enabled: !ui.seriesRunning
                                                text: Number(modelData.heightMin).toFixed(1)
                                                onEditingFinished: Julia.updateSeriesPerturbation(modelData.id, "heightMin", text)
                                            }
                                            Label { text: "Height max" }
                                            TextField {
                                                selectByMouse: true
                                                enabled: !ui.seriesRunning
                                                text: Number(modelData.heightMax).toFixed(1)
                                                onEditingFinished: Julia.updateSeriesPerturbation(modelData.id, "heightMax", text)
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            // -------------------------------------------------- results
            ColumnLayout {
                id: resultsColumn
                property var current: {
                    const results = seriesWindow.results
                    for (let index = 0; index < results.length; ++index) {
                        if (results[index].panel === ui.seriesResultsPanel)
                            return results[index]
                    }
                    return results.length > 0 ? results[0] : null
                }
                Layout.fillWidth: true
                Layout.fillHeight: true
                Layout.margins: 12
                spacing: 8

                RowLayout {
                    Layout.fillWidth: true
                    visible: seriesWindow.results.length > 1
                    spacing: 6

                    Label {
                        text: "Panel"
                        font.bold: true
                    }

                    Repeater {
                        model: seriesWindow.results.length

                        Button {
                            required property int index
                            text: String(index + 1)
                            checkable: true
                            checked: ui.seriesResultsPanel === index + 1
                            highlighted: checked
                            focusPolicy: Qt.NoFocus
                            onClicked: Julia.setSeriesResultsPanel(index + 1)
                        }
                    }

                    Item { Layout.fillWidth: true }
                }

                Rectangle {
                    Layout.fillWidth: true
                    Layout.fillHeight: true
                    visible: resultsColumn.current !== null
                    color: "#ffffff"
                    border.color: "#c4ccd7"
                    border.width: 1
                    radius: 6

                    ColumnLayout {
                        anchors.fill: parent
                        anchors.margins: 10
                        spacing: 8

                        SeriesBarChart {
                            Layout.fillWidth: true
                            Layout.fillHeight: true
                            title: "Head locations"
                            values: resultsColumn.current === null ? [] : resultsColumn.current.locationCounts
                            labels: []
                            xMinimum: resultsColumn.current === null ? NaN : Number(resultsColumn.current.xMin)
                            xMaximum: resultsColumn.current === null ? NaN : Number(resultsColumn.current.xMax)
                            barColor: "#2563eb"
                        }

                        SeriesBarChart {
                            Layout.fillWidth: true
                            Layout.fillHeight: true
                            title: "Number of heads"
                            values: resultsColumn.current === null ? [] : resultsColumn.current.headCounts
                            labels: resultsColumn.current === null ? [] : resultsColumn.current.headLabels
                            barColor: "#7c3aed"
                        }

                        SeriesLineChart {
                            Layout.fillWidth: true
                            Layout.fillHeight: true
                            title: "Patterns"
                            xValues: resultsColumn.current === null ? [] : resultsColumn.current.patternX
                            profiles: resultsColumn.current === null ? [] : resultsColumn.current.patterns
                            converged: resultsColumn.current === null ? [] : resultsColumn.current.patternConverged
                        }
                    }
                }

                Label {
                    Layout.fillWidth: true
                    Layout.fillHeight: true
                    visible: resultsColumn.current === null
                    text: "Results will appear after the first completed realization."
                    color: "#68717d"
                    horizontalAlignment: Text.AlignHCenter
                    verticalAlignment: Text.AlignTop
                }
            }
        }

        Rectangle {
            Layout.fillWidth: true
            visible: ui.message.length > 0
            implicitHeight: seriesMessageText.implicitHeight + 18
            color: "#fff0f0"
            border.color: "#c63f45"

            Label {
                id: seriesMessageText
                anchors.fill: parent
                anchors.margins: 9
                text: ui.message
                color: "#9f252b"
                wrapMode: Text.Wrap
            }
        }
    }
}
