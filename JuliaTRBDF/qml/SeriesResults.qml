import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import QtQuick.Window

Window {
    id: resultWindow

    property string resultsJson: "[]"
    property bool running: false
    property int completedRuns: 0
    property int totalRuns: 0
    property string status: ""
    property var results: {
        try {
            return JSON.parse(resultsJson)
        } catch (error) {
            return []
        }
    }

    width: 1120
    height: 720
    minimumWidth: 720
    minimumHeight: 460
    title: running
           ? "Series results — " + completedRuns + " / " + totalRuns
           : "Series results"
    color: "#eef1f5"

    ColumnLayout {
        anchors.fill: parent
        spacing: 0

        Rectangle {
            Layout.fillWidth: true
            Layout.preferredHeight: 48
            color: "#20252d"

            RowLayout {
                anchors.fill: parent
                anchors.leftMargin: 14
                anchors.rightMargin: 14

                Label {
                    Layout.fillWidth: true
                    text: resultWindow.status
                    color: "white"
                    font.bold: true
                    elide: Text.ElideRight
                }

                Label {
                    text: resultWindow.running
                          ? resultWindow.completedRuns + " / " + resultWindow.totalRuns
                          : ""
                    color: "#bfdbfe"
                }
            }
        }

        ScrollView {
            Layout.fillWidth: true
            Layout.fillHeight: true
            contentWidth: availableWidth
            clip: true

            ColumnLayout {
                width: parent.width
                spacing: 10

                Repeater {
                    model: resultWindow.results

                    Rectangle {
                        required property var modelData
                        Layout.leftMargin: 12
                        Layout.rightMargin: 12
                        Layout.fillWidth: true
                        Layout.preferredHeight: 310
                        color: "#ffffff"
                        border.color: "#c4ccd7"
                        border.width: 1
                        radius: 6

                        ColumnLayout {
                            anchors.fill: parent
                            anchors.margins: 10
                            spacing: 6

                            Label {
                                Layout.fillWidth: true
                                text: "Panel " + modelData.panel
                                color: "#20252d"
                                font.bold: true
                                font.pixelSize: 15
                            }

                            RowLayout {
                                Layout.fillWidth: true
                                Layout.fillHeight: true
                                spacing: 12

                                SeriesBarChart {
                                    Layout.fillWidth: true
                                    Layout.fillHeight: true
                                    title: "Head locations"
                                    subtitle: "local x: " + Number(modelData.xMin).toPrecision(4)
                                              + " — " + Number(modelData.xMax).toPrecision(4)
                                    values: modelData.locationCounts
                                    labels: []
                                    xMinimum: Number(modelData.xMin)
                                    xMaximum: Number(modelData.xMax)
                                    barColor: "#2563eb"
                                }

                                SeriesBarChart {
                                    Layout.fillWidth: true
                                    Layout.fillHeight: true
                                    title: "Number of heads"
                                    subtitle: "per converged realization"
                                    values: modelData.headCounts
                                    labels: modelData.headLabels
                                    barColor: "#7c3aed"
                                }
                            }
                        }
                    }
                }

                Label {
                    Layout.leftMargin: 16
                    Layout.rightMargin: 16
                    Layout.topMargin: 18
                    Layout.fillWidth: true
                    visible: resultWindow.results.length === 0
                    text: "Results will appear after the first completed realization."
                    color: "#68717d"
                    horizontalAlignment: Text.AlignHCenter
                }

                Item { Layout.preferredHeight: 8 }
            }
        }
    }
}
