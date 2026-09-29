import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import QtQuick.Window
import jlqml

ApplicationWindow {
    id: splash
    width: 480
    height: 205
    minimumWidth: width
    maximumWidth: width
    minimumHeight: height
    maximumHeight: height
    visible: true
    title: "Reaction-Diffusion Laboratory"
    // An ordinary window: it can be moved, minimized and covered by others.
    flags: Qt.Window | Qt.CustomizeWindowHint | Qt.WindowTitleHint
           | Qt.WindowMinimizeButtonHint | Qt.WindowCloseButtonHint
    x: Math.round((Screen.width - width) / 2)
    y: Math.round((Screen.height - height) / 2)

    background: Rectangle {
        color: "#f3f5f8"
    }

    // Closing the window cancels the startup, like the Cancel button.
    onClosing: Julia.cancelStartup()

    ColumnLayout {
        anchors.fill: parent
        anchors.margins: 24
        spacing: 12

        Label {
            text: "Loading Reaction-Diffusion Laboratory"
            font.pixelSize: 20
            font.bold: true
            color: "#20252d"
            Layout.fillWidth: true
            // Wider system fonts (e.g. on Linux) shrink instead of clipping.
            fontSizeMode: Text.HorizontalFit
            minimumPixelSize: 14
        }

        Label {
            text: startup.stage
            color: "#4e5a69"
            Layout.fillWidth: true
            elide: Text.ElideRight
        }

        ProgressBar {
            id: bar
            from: 0
            to: 100
            value: startup.percent
            Layout.fillWidth: true

            // Smooths the step to the next stage.
            Behavior on value {
                SmoothedAnimation { velocity: 60 }
            }
        }

        Label {
            text: Math.floor(bar.value) + "%"
            color: "#68717d"
            horizontalAlignment: Text.AlignRight
            Layout.fillWidth: true
        }

        Item { Layout.fillHeight: true }

        Button {
            id: cancelButton
            text: "Cancel"
            Layout.alignment: Qt.AlignRight
            hoverEnabled: true
            onClicked: Julia.cancelStartup()

            // The Basic style has no hover state; darken its own colour instead.
            background: Rectangle {
                implicitWidth: 100
                implicitHeight: 40
                color: cancelButton.down ? Qt.darker(cancelButton.palette.button, 1.25)
                     : cancelButton.hovered ? Qt.darker(cancelButton.palette.button, 1.12)
                     : cancelButton.palette.button
                border.color: cancelButton.palette.highlight
                border.width: cancelButton.visualFocus ? 2 : 0

                Behavior on color {
                    ColorAnimation { duration: 80 }
                }
            }
        }
    }
}
