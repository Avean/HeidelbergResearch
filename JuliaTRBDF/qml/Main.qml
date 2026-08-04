import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import QtQml.Models
import jlqml
import Makie

ApplicationWindow {
    id: window

    visible: true
    width: 1500
    height: 900
    minimumWidth: 980
    minimumHeight: 650
    title: "Reaction-Diffusion Laboratory"
    color: "#eef1f5"

    property var modelCatalog: JSON.parse(ui.modelCatalogJson)
    property var variables: JSON.parse(ui.variablesJson)
    property var equationImages: JSON.parse(ui.equationImagesJson)
    property bool textEditorFocused: false
    property int selectedFamilyIndex: 0
    property int activeFamilyIndex: findFamilyIndex(ui.activeModelKey)
    property string bottomPanel: ""

    function findFamilyIndex(modelKey) {
        for (let familyIndex = 0; familyIndex < modelCatalog.length; ++familyIndex) {
            const models = modelCatalog[familyIndex].models
            for (let modelIndex = 0; modelIndex < models.length; ++modelIndex) {
                if (models[modelIndex].key === modelKey)
                    return familyIndex
            }
        }
        return 0
    }

    function selectedFamilyModels() {
        if (selectedFamilyIndex < 0 || selectedFamilyIndex >= modelCatalog.length)
            return []
        return modelCatalog[selectedFamilyIndex].models
    }

    function activeModelIndexInSelectedFamily() {
        const models = selectedFamilyModels()
        for (let index = 0; index < models.length; ++index) {
            if (models[index].key === ui.activeModelKey)
                return index
        }
        return -1
    }

    function openModelDrawer() {
        controlDrawer.close()
        modelDrawer.open()
    }

    function openControlDrawer() {
        modelDrawer.close()
        controlDrawer.open()
    }

    function toggleBottomPanel(panelName) {
        bottomPanel = bottomPanel === panelName ? "" : panelName
    }

    onActiveFamilyIndexChanged: selectedFamilyIndex = activeFamilyIndex

    Component.onCompleted: selectedFamilyIndex = activeFamilyIndex

    palette.window: "#eef1f5"
    palette.windowText: "#20252d"
    palette.button: "#f5f7fa"
    palette.buttonText: "#20252d"
    palette.highlight: "#3b82f6"
    palette.highlightedText: "white"

    header: ToolBar {
        id: topBar
        height: 52

        background: Rectangle {
            color: "#20252d"
        }

        RowLayout {
            anchors.fill: parent
            anchors.leftMargin: 8
            anchors.rightMargin: 8
            spacing: 8

            ToolButton {
                text: "Models: " + ui.modelName
                enabled: !ui.graphicsBusy
                palette.buttonText: "white"
                onClicked: modelDrawer.opened ? modelDrawer.close() : window.openModelDrawer()
            }

            Item {
                Layout.fillWidth: true
            }

            Rectangle {
                id: runningButton
                Layout.preferredWidth: 96
                Layout.preferredHeight: 34
                radius: 5
                color: ui.running ? "#26945b" : "#c63f45"

                Text {
                    anchors.centerIn: parent
                    text: ui.running ? "Running" : "Stopped"
                    color: "white"
                    font.bold: true
                }

                MouseArea {
                    anchors.fill: parent
                    enabled: !ui.graphicsBusy
                    cursorShape: Qt.PointingHandCursor
                    onClicked: Julia.toggleRunning()
                }
            }

            ToolButton {
                text: controlDrawer.opened ? "Close controls" : "Controls"
                palette.buttonText: "white"
                onClicked: controlDrawer.opened ? controlDrawer.close() : window.openControlDrawer()
            }
        }
    }

    footer: ToolBar {
        id: bottomBar
        height: 48

        background: Rectangle {
            color: "#20252d"
        }

        RowLayout {
            anchors.fill: parent
            anchors.leftMargin: 12
            anchors.rightMargin: 12
            spacing: 10

            Label {
                text: "Domain rescale"
                color: "white"
                font.bold: true
            }

            Slider {
                Layout.preferredWidth: Math.min(330, window.width * 0.27)
                enabled: !ui.graphicsBusy
                from: 0
                to: 3
                stepSize: 0.2
                value: 2 * Math.log(Number(ui.domainLength)) / Math.LN10
                onMoved: Julia.setDomainExponent(value)
            }

            Label {
                Layout.preferredWidth: 72
                text: Number(ui.domainLength).toPrecision(3)
                color: "white"
                font.family: "Consolas"
            }

            Item {
                Layout.fillWidth: true
            }

            ToolButton {
                text: "Perturbations"
                palette.buttonText: window.bottomPanel === "perturbations" ? "#93c5fd" : "white"
                onClicked: window.toggleBottomPanel("perturbations")
            }

            ToolButton {
                text: "Split / Merge"
                palette.buttonText: window.bottomPanel === "partition" ? "#93c5fd" : "white"
                onClicked: window.toggleBottomPanel("partition")
            }
        }
    }

    MakieArea {
        id: plotArea
        anchors.fill: parent
        scene: plot
    }

    Rectangle {
        id: bottomDrawer
        z: 35
        anchors.left: parent.left
        anchors.right: parent.right
        anchors.bottom: parent.bottom
        height: window.bottomPanel === "" ? 0
                : window.bottomPanel === "perturbations" ? 112 : 150
        color: "#f3f5f8"
        border.color: "#929ca9"
        border.width: height > 0 ? 1 : 0
        clip: true
        visible: height > 0

        Behavior on height {
            NumberAnimation {
                duration: 180
                easing.type: Easing.OutCubic
            }
        }

        StackLayout {
            anchors.fill: parent
            currentIndex: window.bottomPanel === "partition" ? 1 : 0

            Item {
                RowLayout {
                    anchors.fill: parent
                    anchors.leftMargin: 16
                    anchors.rightMargin: 10
                    anchors.topMargin: 10
                    anchors.bottomMargin: 10
                    spacing: 10

                    Label {
                        text: "Perturbations"
                        font.bold: true
                        font.pixelSize: 16
                    }

                    Button {
                        text: (ui.randomMode ? "Random" : "Constant") + "  [Z]"
                        highlighted: ui.randomMode
                        enabled: !ui.graphicsBusy
                        onClicked: Julia.toggleRandomMode()
                    }

                    Button {
                        text: (ui.absoluteMode ? "Absolute" : "Relative") + "  [X]"
                        highlighted: ui.absoluteMode
                        enabled: !ui.graphicsBusy
                        onClicked: Julia.toggleAbsoluteMode()
                    }

                    Label {
                        text: "Width"
                    }

                    TextField {
                        id: perturbationWidthField
                        Layout.preferredWidth: 82
                        selectByMouse: true
                        validator: DoubleValidator {
                            bottom: 0.0000000001
                            top: 1.0
                            notation: DoubleValidator.StandardNotation
                        }
                        onActiveFocusChanged: window.textEditorFocused = activeFocus
                        onEditingFinished: Julia.setPerturbationWidth(text)

                        Binding on text {
                            value: Number(ui.perturbationWidth).toFixed(2)
                            when: !perturbationWidthField.activeFocus
                            restoreMode: Binding.RestoreBindingOrValue
                        }
                    }

                    Label {
                        visible: ui.absoluteMode
                        text: "Height"
                    }

                    TextField {
                        id: perturbationHeightField
                        visible: ui.absoluteMode
                        Layout.preferredWidth: 92
                        selectByMouse: true
                        validator: DoubleValidator {
                            notation: DoubleValidator.ScientificNotation
                        }
                        onActiveFocusChanged: window.textEditorFocused = activeFocus
                        onEditingFinished: Julia.setPerturbationHeight(text)

                        Binding on text {
                            value: Number(ui.perturbationHeight).toString()
                            when: !perturbationHeightField.activeFocus
                            restoreMode: Binding.RestoreBindingOrValue
                        }
                    }

                    Label {
                        Layout.fillWidth: true
                        text: ui.absoluteMode
                              ? "Scroll: width   |   Ctrl + scroll: height"
                              : "Scroll: width   |   Ctrl + scroll: relative preview scale"
                        color: "#68717d"
                        horizontalAlignment: Text.AlignHCenter
                        elide: Text.ElideRight
                    }

                    ToolButton {
                        text: "Close"
                        onClicked: window.bottomPanel = ""
                    }
                }
            }

            Item {
                ColumnLayout {
                    anchors.fill: parent
                    anchors.leftMargin: 12
                    anchors.rightMargin: 10
                    anchors.topMargin: 8
                    anchors.bottomMargin: 8
                    spacing: 6

                    RowLayout {
                        Layout.fillWidth: true
                        spacing: 8

                        Label {
                            text: "Split / Merge"
                            font.bold: true
                            font.pixelSize: 16
                        }

                        Button {
                            text: "Left"
                            enabled: !ui.graphicsBusy && ui.selectedSegment > 1
                            onClicked: Julia.changeSelectedSegment(-1)
                        }

                        Label {
                            Layout.preferredWidth: 100
                            horizontalAlignment: Text.AlignHCenter
                            text: "Panel " + ui.selectedSegment + " / " + ui.segmentCount
                        }

                        Button {
                            text: "Right"
                            enabled: !ui.graphicsBusy && ui.selectedSegment < ui.segmentCount
                            onClicked: Julia.changeSelectedSegment(1)
                        }

                        Item {
                            Layout.fillWidth: true
                        }

                        Rectangle {
                            Layout.preferredWidth: 136
                            Layout.preferredHeight: 34
                            visible: ui.segmentCount > 1
                            radius: 5
                            color: ui.synchronizationStatus === "Synchronized"
                                   ? "#26945b"
                                   : ui.synchronizationStatus === "Synchronizing..."
                                     ? "#d79a22" : "#e06a2f"

                            Text {
                                anchors.centerIn: parent
                                text: ui.synchronizationStatus
                                color: "white"
                                font.bold: true
                            }

                            MouseArea {
                                anchors.fill: parent
                                enabled: !ui.graphicsBusy
                                cursorShape: Qt.PointingHandCursor
                                onClicked: Julia.synchronizeDomains()
                            }
                        }

                        ToolButton {
                            text: "Close"
                            onClicked: window.bottomPanel = ""
                        }
                    }

                    RowLayout {
                        Layout.fillWidth: true
                        spacing: 8

                        Label {
                            text: "Split point: " + ui.splitIndex
                        }

                        Slider {
                            Layout.preferredWidth: Math.min(360, window.width * 0.28)
                            enabled: !ui.graphicsBusy
                            from: 2
                            to: Math.max(2, ui.splitMaximum)
                            stepSize: 1
                            value: ui.splitIndex
                            onMoved: Julia.setSplitIndex(Math.round(value))
                        }

                        Button {
                            enabled: !ui.graphicsBusy
                            text: ui.graphicsBusy ? "Updating..." : "Split selected panel"
                            onClicked: Julia.splitSelectedSegment()
                        }

                        Label {
                            visible: ui.segmentCount > 1
                            text: "Merge:"
                            font.bold: true
                        }

                        ScrollView {
                            Layout.fillWidth: true
                            Layout.preferredHeight: 42
                            visible: ui.segmentCount > 1
                            contentHeight: availableHeight
                            ScrollBar.vertical.policy: ScrollBar.AlwaysOff
                            ScrollBar.horizontal.policy: ScrollBar.AsNeeded
                            clip: true

                            Row {
                                spacing: 6

                                Repeater {
                                    model: Math.max(0, ui.segmentCount - 1)

                                    Button {
                                        required property int index
                                        enabled: !ui.graphicsBusy
                                        text: (index + 1) + " | " + (index + 2)
                                        onClicked: Julia.mergeBoundary(index + 1)
                                    }
                                }
                            }
                        }

                        Label {
                            Layout.fillWidth: true
                            visible: ui.segmentCount <= 1
                            text: "No divided domains to merge or synchronize."
                            color: "#68717d"
                            horizontalAlignment: Text.AlignHCenter
                        }
                    }
                }
            }
        }
    }

    Rectangle {
        id: leftEdgeHotspot
        z: 20
        width: 9
        color: "transparent"
        anchors.top: parent.top
        anchors.bottom: parent.bottom
        anchors.left: parent.left
        visible: !modelDrawer.opened

        HoverHandler {
            id: leftEdgeHover
            acceptedDevices: PointerDevice.Mouse | PointerDevice.TouchPad
            onHoveredChanged: hovered ? modelOpenDelay.restart() : modelOpenDelay.stop()
        }
    }

    Timer {
        id: modelOpenDelay
        interval: 350
        repeat: false
        onTriggered: {
            if (leftEdgeHover.hovered && !modelDrawer.opened)
                window.openModelDrawer()
        }
    }

    Timer {
        id: modelCloseDelay
        interval: 500
        repeat: false
        onTriggered: {
            if (!modelDrawerHover.hovered && !modelDrawer.pinned)
                modelDrawer.close()
        }
    }

    Drawer {
        id: modelDrawer
        property bool pinned: false

        edge: Qt.LeftEdge
        width: Math.min(530, window.width * 0.46)
        height: window.height - topBar.height - bottomBar.height
        y: topBar.height
        modal: false
        dim: false
        interactive: true
        closePolicy: Popup.NoAutoClose

        background: Rectangle {
            color: "#f3f5f8"
            border.color: "#aab3bf"
            border.width: 1
        }

        contentItem: Rectangle {
            color: "#f3f5f8"

            HoverHandler {
                id: modelDrawerHover
                acceptedDevices: PointerDevice.Mouse | PointerDevice.TouchPad
                onHoveredChanged: {
                    if (hovered)
                        modelCloseDelay.stop()
                    else if (!modelDrawer.pinned)
                        modelCloseDelay.restart()
                }
            }

            ColumnLayout {
                anchors.fill: parent
                spacing: 0

                Rectangle {
                    Layout.fillWidth: true
                    Layout.preferredHeight: 48
                    color: "#2b313b"

                    RowLayout {
                        anchors.fill: parent
                        anchors.leftMargin: 14
                        anchors.rightMargin: 8

                        Label {
                            text: "Models and equations"
                            color: "white"
                            font.bold: true
                            font.pixelSize: 16
                            Layout.fillWidth: true
                        }

                        ToolButton {
                            text: modelDrawer.pinned ? "Unpin" : "Pin"
                            palette.buttonText: "white"
                            onClicked: modelDrawer.pinned = !modelDrawer.pinned
                        }

                        ToolButton {
                            text: "Close"
                            palette.buttonText: "white"
                            onClicked: modelDrawer.close()
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
                        spacing: 9

                        Item {
                            Layout.preferredHeight: 2
                        }

                        ControlSection {
                            title: "Model selection"
                            Layout.leftMargin: 9
                            Layout.rightMargin: 9

                            Label {
                                Layout.fillWidth: true
                                text: "Model family"
                            }

                            ComboBox {
                                id: familyCombo
                                Layout.fillWidth: true
                                enabled: !ui.graphicsBusy
                                model: window.modelCatalog.map(function(item) { return item.family })
                                currentIndex: window.selectedFamilyIndex
                                onActivated: window.selectedFamilyIndex = currentIndex
                            }

                            Label {
                                Layout.fillWidth: true
                                text: "Model"
                            }

                            ComboBox {
                                id: modelCombo
                                property var entries: window.selectedFamilyModels()
                                Layout.fillWidth: true
                                enabled: !ui.graphicsBusy
                                model: entries.map(function(item) { return item.label })
                                currentIndex: window.activeModelIndexInSelectedFamily()
                                displayText: currentIndex >= 0 ? currentText : "Select model"
                                onActivated: Julia.selectModel(entries[currentIndex].key)
                            }
                        }

                        ControlSection {
                            title: "Equations"
                            Layout.leftMargin: 9
                            Layout.rightMargin: 9

                            Repeater {
                                model: window.equationImages

                                Image {
                                    required property var modelData
                                    Layout.fillWidth: true
                                    Layout.preferredHeight: implicitWidth > 0
                                                            ? Math.max(42, Math.min(190, width * implicitHeight / implicitWidth))
                                                            : 60
                                    source: modelData
                                    fillMode: Image.PreserveAspectFit
                                    horizontalAlignment: Image.AlignLeft
                                    asynchronous: false
                                    cache: true
                                }
                            }

                            Label {
                                Layout.fillWidth: true
                                visible: window.equationImages.length === 0
                                text: "No equations specified for this model."
                                color: "#68717d"
                            }
                        }

                        ControlSection {
                            title: "Boundary conditions"
                            Layout.leftMargin: 9
                            Layout.rightMargin: 9

                            RowLayout {
                                Layout.fillWidth: true

                                Button {
                                    Layout.fillWidth: true
                                    text: "Neumann"
                                    highlighted: ui.boundaryName === text
                                    enabled: !ui.graphicsBusy
                                    onClicked: Julia.selectBoundaryCondition(text)
                                }

                                Button {
                                    Layout.fillWidth: true
                                    text: "Periodic"
                                    highlighted: ui.boundaryName === text
                                    enabled: !ui.graphicsBusy
                                    onClicked: Julia.selectBoundaryCondition(text)
                                }
                            }
                        }

                        Item {
                            Layout.preferredHeight: 8
                        }
                    }
                }
            }
        }
    }

    Rectangle {
        id: rightEdgeHotspot
        z: 20
        width: 9
        color: "transparent"
        anchors.top: parent.top
        anchors.bottom: parent.bottom
        anchors.right: parent.right
        visible: !controlDrawer.opened

        HoverHandler {
            id: rightEdgeHover
            acceptedDevices: PointerDevice.Mouse | PointerDevice.TouchPad
            onHoveredChanged: hovered ? controlOpenDelay.restart() : controlOpenDelay.stop()
        }
    }

    Timer {
        id: controlOpenDelay
        interval: 350
        repeat: false
        onTriggered: {
            if (rightEdgeHover.hovered && !controlDrawer.opened)
                window.openControlDrawer()
        }
    }

    Timer {
        id: controlCloseDelay
        interval: 500
        repeat: false
        onTriggered: {
            if (!controlDrawerHover.hovered && !controlDrawer.pinned)
                controlDrawer.close()
        }
    }

    Drawer {
        id: controlDrawer
        property bool pinned: false

        edge: Qt.RightEdge
        width: Math.min(440, window.width * 0.42)
        height: window.height - topBar.height - bottomBar.height
        y: topBar.height
        modal: false
        dim: false
        interactive: true
        closePolicy: Popup.NoAutoClose

        background: Rectangle {
            color: "#f3f5f8"
            border.color: "#aab3bf"
            border.width: 1
        }

        contentItem: Rectangle {
            color: "#f3f5f8"

            HoverHandler {
                id: controlDrawerHover
                acceptedDevices: PointerDevice.Mouse | PointerDevice.TouchPad
                onHoveredChanged: {
                    if (hovered)
                        controlCloseDelay.stop()
                    else if (!controlDrawer.pinned)
                        controlCloseDelay.restart()
                }
            }

            ColumnLayout {
                anchors.fill: parent
                spacing: 0

                Rectangle {
                    Layout.fillWidth: true
                    Layout.preferredHeight: 48
                    color: "#2b313b"

                    RowLayout {
                        anchors.fill: parent
                        anchors.leftMargin: 14
                        anchors.rightMargin: 8

                        Label {
                            text: "Simulation controls"
                            color: "white"
                            font.bold: true
                            font.pixelSize: 16
                            Layout.fillWidth: true
                        }

                        ToolButton {
                            text: controlDrawer.pinned ? "Unpin" : "Pin"
                            palette.buttonText: "white"
                            onClicked: controlDrawer.pinned = !controlDrawer.pinned
                        }

                        ToolButton {
                            text: "Close"
                            palette.buttonText: "white"
                            onClicked: controlDrawer.close()
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
                        spacing: 9

                        Item {
                            Layout.preferredHeight: 2
                        }

                        ControlSection {
                            title: "Time and reset"
                            Layout.leftMargin: 9
                            Layout.rightMargin: 9

                            Label {
                                Layout.fillWidth: true
                                text: "Maximum time step: " + Number(ui.dtmax).toExponential(1)
                            }

                            Slider {
                                Layout.fillWidth: true
                                enabled: !ui.graphicsBusy
                                from: -5
                                to: 5
                                stepSize: 1
                                value: Math.log(Number(ui.dtmax)) / Math.LN10
                                onMoved: Julia.setDtExponent(Math.round(value))
                            }

                            Button {
                                Layout.fillWidth: true
                                enabled: !ui.graphicsBusy
                                text: ui.graphicsBusy ? "Updating plots..." : "Reset initial state"
                                onClicked: Julia.resetSimulation()
                            }
                        }

                        ControlSection {
                            title: "Constant initial condition"
                            Layout.leftMargin: 9
                            Layout.rightMargin: 9

                            Repeater {
                                model: window.variables

                                RowLayout {
                                    required property int index
                                    required property var modelData
                                    Layout.fillWidth: true

                                    Label {
                                        text: modelData
                                        Layout.preferredWidth: 45
                                    }

                                    TextField {
                                        id: constantValue
                                        Layout.fillWidth: true
                                        text: "0.0"
                                        selectByMouse: true
                                        validator: DoubleValidator {}
                                        onActiveFocusChanged: window.textEditorFocused = activeFocus
                                    }

                                    Button {
                                        text: "Apply"
                                        enabled: !ui.graphicsBusy
                                        onClicked: Julia.applyConstantInitialCondition(index, constantValue.text)
                                    }
                                }
                            }
                        }

                        Item {
                            Layout.preferredHeight: 8
                        }
                    }
                }
            }
        }
    }

    Rectangle {
        id: messageBanner
        z: 60
        visible: ui.message.length > 0
        anchors.horizontalCenter: parent.horizontalCenter
        anchors.bottom: bottomDrawer.height > 0 ? bottomDrawer.top : parent.bottom
        width: Math.min(parent.width - 40, 760)
        implicitHeight: messageText.implicitHeight + 18
        color: "#fff0f0"
        border.color: "#c63f45"
        radius: 5

        Label {
            id: messageText
            anchors.fill: parent
            anchors.margins: 9
            text: ui.message
            color: "#9f252b"
            wrapMode: Text.Wrap
        }
    }

    Shortcut {
        sequence: "Z"
        context: Qt.WindowShortcut
        enabled: window.active && !window.textEditorFocused && !ui.graphicsBusy
        autoRepeat: false
        onActivated: Julia.toggleRandomMode()
    }

    Shortcut {
        sequence: "X"
        context: Qt.WindowShortcut
        enabled: window.active && !window.textEditorFocused && !ui.graphicsBusy
        autoRepeat: false
        onActivated: Julia.toggleAbsoluteMode()
    }

    Timer {
        interval: 33
        running: true
        repeat: true
        onTriggered: {
            Julia.refreshUI()
            plotArea.update()
        }
    }

    Timer {
        interval: Math.max(1, ui.autoCloseMs)
        running: ui.autoCloseMs > 0
        repeat: false
        onTriggered: window.close()
    }

    onClosing: function(close) {
        Julia.requestClose()
    }
}
