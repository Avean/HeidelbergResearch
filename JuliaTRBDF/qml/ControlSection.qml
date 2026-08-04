import QtQuick
import QtQuick.Controls
import QtQuick.Layouts

Rectangle {
    id: root

    property string title: "Section"
    property bool expanded: true
    default property alias sectionContent: contentColumn.data

    Layout.fillWidth: true
    implicitHeight: sectionLayout.implicitHeight
    color: "#f8fafc"
    border.color: "#d7dde5"
    border.width: 1
    radius: 6

    ColumnLayout {
        id: sectionLayout
        anchors.left: parent.left
        anchors.right: parent.right
        anchors.top: parent.top
        spacing: 0

        ToolButton {
            Layout.fillWidth: true
            text: (root.expanded ? "▼  " : "▶  ") + root.title
            font.bold: true
            leftPadding: 12
            rightPadding: 12
            contentItem: Text {
                text: parent.text
                color: "#20252d"
                font: parent.font
                horizontalAlignment: Text.AlignLeft
                verticalAlignment: Text.AlignVCenter
            }
            background: Rectangle {
                color: parent.hovered ? "#e9eef5" : "transparent"
                radius: 6
            }
            onClicked: root.expanded = !root.expanded
        }

        ColumnLayout {
            id: contentColumn
            Layout.fillWidth: true
            Layout.leftMargin: 12
            Layout.rightMargin: 12
            Layout.bottomMargin: root.expanded ? 12 : 0
            spacing: 8
            visible: root.expanded
        }
    }
}
