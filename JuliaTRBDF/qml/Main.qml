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

    ButtonGroup {
        id: splitSegmentButtonGroup
        exclusive: true
    }

    ButtonGroup {
        id: steadySegmentButtonGroup
        exclusive: true
    }

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

    function dtMatchesExponent(exponent) {
        const dt = Number(ui.dtmax)
        return dt > 0
                && Math.abs(Math.log(dt) / Math.LN10 - exponent) < 0.0001
    }

    function openModelDrawer() {
        controlDrawer.close()
        bottomDrawer.close()
        modelDrawer.open()
    }

    function openControlDrawer() {
        modelDrawer.close()
        bottomDrawer.close()
        controlDrawer.open()
    }

    function toggleBottomPanel(panelName) {
        modelDrawer.close()
        controlDrawer.close()

        if (bottomDrawer.opened && bottomPanel === panelName) {
            closeBottomPanel()
        } else {
            bottomPanel = panelName
            if (!bottomDrawer.opened)
                bottomDrawer.open()
        }
    }

    function closeBottomPanel() {
        perturbationWidthField.focus = false
        perturbationHeightField.focus = false
        textEditorFocused = false
        bottomDrawer.close()
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
            spacing: 6

            ToolButton {
                text: "Models: " + ui.modelName
                Layout.maximumWidth: Math.min(390, window.width * 0.30)
                enabled: !ui.graphicsBusy
                palette.buttonText: "white"
                onClicked: modelDrawer.opened ? modelDrawer.close() : window.openModelDrawer()
            }

            ToolButton {
                id: stateButton
                text: ui.checkpointAvailable ? "State  •" : "State"
                enabled: !ui.graphicsBusy
                palette.buttonText: "white"
                onClicked: stateMenu.open()

                Menu {
                    id: stateMenu
                    y: stateButton.height

                    MenuItem {
                        text: "Save current state"
                        enabled: !ui.graphicsBusy
                        onTriggered: Julia.saveCurrentState()
                    }

                    MenuItem {
                        text: "Restore saved state"
                        enabled: ui.checkpointAvailable && !ui.graphicsBusy
                        onTriggered: Julia.restoreSavedState()
                    }
                }
            }

            Item {
                Layout.fillWidth: true
            }

            Rectangle {
                id: runningButton
                Layout.preferredWidth: 92
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

            Item {
                Layout.preferredWidth: 5
            }

            Label {
                text: "Speed"
                color: "white"
                font.bold: true
            }

            Repeater {
                model: [
                    { key: "1", exponent: -3, description: "Slow: maximum dt = 1e-3" },
                    { key: "2", exponent: -1, description: "Medium: maximum dt = 1e-1" },
                    { key: "3", exponent: 1, description: "Fast: maximum dt = 1e1" },
                    { key: "4", exponent: 5, description: "Very fast: maximum dt = 1e5" }
                ]

                Rectangle {
                    id: speedPreset
                    required property var modelData
                    Layout.preferredWidth: 30
                    Layout.preferredHeight: 30
                    property bool active: window.dtMatchesExponent(modelData.exponent)
                    radius: 4
                    color: active ? "#3b82f6" : "#46515f"
                    border.color: active ? "#93c5fd" : "#697586"
                    opacity: ui.graphicsBusy ? 0.55 : 1.0

                    Text {
                        anchors.centerIn: parent
                        text: speedPreset.modelData.key
                        color: "white"
                        font.bold: true
                    }

                    MouseArea {
                        id: speedPresetMouse
                        anchors.fill: parent
                        enabled: !ui.graphicsBusy
                        hoverEnabled: true
                        cursorShape: Qt.PointingHandCursor
                        onClicked: Julia.setDtExponent(speedPreset.modelData.exponent)
                    }

                    ToolTip.visible: speedPresetMouse.containsMouse
                    ToolTip.text: modelData.description + "  [" + modelData.key + "]"
                }
            }

            Label {
                text: window.width < 1150 ? "Max dt" : "Maximum dt"
                color: "white"
                font.bold: true
            }

            Slider {
                Layout.preferredWidth: Math.min(
                    210,
                    Math.max(110, window.width * 0.13)
                )
                enabled: !ui.graphicsBusy
                from: -5
                to: 5
                stepSize: 1
                value: Math.log(Number(ui.dtmax)) / Math.LN10
                onMoved: Julia.setDtExponent(Math.round(value))
            }

            Label {
                Layout.preferredWidth: 62
                text: Number(ui.dtmax).toExponential(1)
                color: "white"
                font.family: "Consolas"
            }

            Rectangle {
                Layout.preferredWidth: 68
                Layout.preferredHeight: 34
                radius: 5
                color: "#596575"

                Text {
                    anchors.centerIn: parent
                    text: ui.graphicsBusy ? "Wait..." : "Reset"
                    color: "white"
                    font.bold: true
                }

                MouseArea {
                    anchors.fill: parent
                    enabled: !ui.graphicsBusy
                    cursorShape: Qt.PointingHandCursor
                    onClicked: Julia.resetSimulation()
                }
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
            spacing: 8

            Item {
                Layout.fillWidth: true
                Layout.preferredWidth: 1
                Layout.minimumWidth: 230

                RowLayout {
                    anchors.left: parent.left
                    anchors.verticalCenter: parent.verticalCenter
                    spacing: 6

                    ToolButton {
                        Layout.alignment: Qt.AlignVCenter
                        text: "Split / Merge"
                        palette.buttonText: window.bottomPanel === "partition" ? "#93c5fd" : "white"
                        onClicked: window.toggleBottomPanel("partition")
                    }

                    ToolButton {
                        Layout.alignment: Qt.AlignVCenter
                        text: "Perturbations"
                        palette.buttonText: window.bottomPanel === "perturbations" ? "#93c5fd" : "white"
                        onClicked: window.toggleBottomPanel("perturbations")
                    }
                }
            }

            RowLayout {
                Layout.alignment: Qt.AlignHCenter
                spacing: 8

                Label {
                    text: "Domain rescale"
                    color: "white"
                    font.bold: true
                }

                Slider {
                    Layout.preferredWidth: Math.min(300, window.width * 0.24)
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
            }

            Item {
                Layout.fillWidth: true
                Layout.preferredWidth: 1
                Layout.minimumWidth: 230

                ToolButton {
                    anchors.right: parent.right
                    anchors.verticalCenter: parent.verticalCenter
                    text: controlDrawer.opened ? "Close steady state" : "Set steady state"
                    palette.buttonText: "white"
                    onClicked: controlDrawer.opened ? controlDrawer.close() : window.openControlDrawer()
                }
            }
        }
    }

    MakieArea {
        id: plotArea
        anchors.fill: parent
        scene: plot
    }

    Popup {
        id: bottomDrawer
        x: 0
        y: parent.height - height
        width: parent.width
        height: window.bottomPanel === "perturbations" ? 108 : 174
        modal: true
        dim: false
        focus: true
        padding: 0
        closePolicy: Popup.CloseOnPressOutside | Popup.CloseOnEscape

        enter: Transition {
            NumberAnimation {
                property: "y"
                from: bottomDrawer.parent.height
                to: bottomDrawer.parent.height - bottomDrawer.height
                duration: 180
                easing.type: Easing.OutCubic
            }
        }

        exit: Transition {
            NumberAnimation {
                property: "y"
                from: bottomDrawer.parent.height - bottomDrawer.height
                to: bottomDrawer.parent.height
                duration: 150
                easing.type: Easing.InCubic
            }
        }

        onClosed: {
            perturbationWidthField.focus = false
            perturbationHeightField.focus = false
            window.textEditorFocused = false
            window.bottomPanel = ""
        }

        background: Rectangle {
            color: "#f3f5f8"
            border.color: "#929ca9"
            border.width: 1
        }

        Behavior on height {
            NumberAnimation {
                duration: 180
                easing.type: Easing.OutCubic
            }
        }

        contentItem: StackLayout {
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
                        onClicked: window.closeBottomPanel()
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

                        Label {
                            text: "Panel"
                            font.bold: true
                        }

                        ScrollView {
                            Layout.preferredWidth: Math.min(
                                330,
                                Math.max(48, ui.segmentCount * 50)
                            )
                            Layout.preferredHeight: 38
                            contentHeight: availableHeight
                            ScrollBar.vertical.policy: ScrollBar.AlwaysOff
                            ScrollBar.horizontal.policy: ScrollBar.AsNeeded
                            clip: true

                            Row {
                                spacing: 6

                                Repeater {
                                    model: ui.segmentCount

                                    Button {
                                        required property int index
                                        width: 44
                                        height: 32
                                        text: String(index + 1)
                                        checkable: true
                                        checked: ui.selectedSegment === index + 1
                                        highlighted: checked
                                        focusPolicy: Qt.NoFocus
                                        enabled: !ui.graphicsBusy
                                        ButtonGroup.group: splitSegmentButtonGroup
                                        onClicked: Julia.selectSplitSegment(index + 1)
                                    }
                                }
                            }
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
                            onClicked: window.closeBottomPanel()
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

                        Button {
                            enabled: ui.segmentCount > 1 && !ui.graphicsBusy
                            text: "Delete selected panel"
                            onClicked: Julia.deleteSelectedSegment()
                        }

                        Item {
                            Layout.fillWidth: true
                        }

                        Label {
                            visible: ui.segmentCount <= 1
                            text: "Split the domain to enable Merge, Swap and Delete."
                            color: "#68717d"
                        }
                    }

                    RowLayout {
                        Layout.fillWidth: true
                        visible: ui.segmentCount > 1
                        spacing: 8

                        Label {
                            text: "Merge:"
                            font.bold: true
                        }

                        ScrollView {
                            Layout.fillWidth: true
                            Layout.preferredHeight: 42
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
                            text: "Swap:"
                            font.bold: true
                        }

                        ScrollView {
                            Layout.fillWidth: true
                            Layout.preferredHeight: 42
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
                                        text: (index + 1) + " ↔ " + (index + 2)
                                        onClicked: Julia.swapBoundary(index + 1)
                                    }
                                }
                            }
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

    Drawer {
        id: modelDrawer

        edge: Qt.LeftEdge
        width: Math.min(
            window.width * 0.90,
            Math.max(410, Number(ui.equationPreferredWidth) + 58)
        )
        height: window.height - topBar.height - bottomBar.height
        y: topBar.height
        modal: true
        dim: false
        interactive: true
        closePolicy: Popup.CloseOnPressOutside | Popup.CloseOnEscape

        Behavior on width {
            NumberAnimation {
                duration: 180
                easing.type: Easing.OutCubic
            }
        }

        background: Rectangle {
            color: "#f3f5f8"
            border.color: "#aab3bf"
            border.width: 1
        }

        contentItem: Rectangle {
            color: "#f3f5f8"

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
                                popup.height: Math.min(
                                    popup.implicitContentHeight + popup.topPadding + popup.bottomPadding,
                                    320
                                )
                                onActivated: {
                                    window.selectedFamilyIndex = currentIndex
                                    const entries = window.selectedFamilyModels()
                                    if (entries.length > 0)
                                        Julia.selectModel(entries[0].key)
                                }
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
                                popup.height: Math.min(
                                    popup.implicitContentHeight + popup.topPadding + popup.bottomPadding,
                                    360
                                )
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
                                                            ? Math.max(90, Math.min(620, width * implicitHeight / implicitWidth))
                                                            : 90
                                    source: modelData
                                    sourceSize.width: Math.max(
                                        2400,
                                        Math.ceil(width * window.screen.devicePixelRatio * 2.5)
                                    )
                                    fillMode: Image.PreserveAspectFit
                                    horizontalAlignment: Image.AlignLeft
                                    asynchronous: false
                                    cache: true
                                    smooth: true
                                    mipmap: true
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

    Popup {
        id: controlDrawer
        x: 0
        y: parent.height - height
        width: parent.width
        height: 132
        modal: true
        dim: false
        focus: true
        padding: 0
        closePolicy: Popup.CloseOnPressOutside | Popup.CloseOnEscape

        enter: Transition {
            NumberAnimation {
                property: "y"
                from: controlDrawer.parent.height
                to: controlDrawer.parent.height - controlDrawer.height
                duration: 180
                easing.type: Easing.OutCubic
            }
        }

        exit: Transition {
            NumberAnimation {
                property: "y"
                from: controlDrawer.parent.height - controlDrawer.height
                to: controlDrawer.parent.height
                duration: 150
                easing.type: Easing.InCubic
            }
        }

        onClosed: window.textEditorFocused = false

        background: Rectangle {
            color: "#f3f5f8"
            border.color: "#aab3bf"
            border.width: 1
        }

        contentItem: ColumnLayout {
            anchors.fill: parent
            spacing: 0

            Rectangle {
                Layout.fillWidth: true
                Layout.preferredHeight: 46
                color: "#2b313b"

                RowLayout {
                    anchors.fill: parent
                    anchors.leftMargin: 14
                    anchors.rightMargin: 8
                    spacing: 9

                    Label {
                        text: "Set steady state"
                        color: "white"
                        font.bold: true
                        font.pixelSize: 16
                    }

                    Label {
                        text: "Target panel"
                        color: "#cbd3dc"
                        font.bold: true
                    }

                    ScrollView {
                        Layout.preferredWidth: Math.min(
                            330,
                            Math.max(48, ui.segmentCount * 50)
                        )
                        Layout.preferredHeight: 38
                        contentHeight: availableHeight
                        ScrollBar.vertical.policy: ScrollBar.AlwaysOff
                        ScrollBar.horizontal.policy: ScrollBar.AsNeeded
                        clip: true

                        Row {
                            spacing: 6

                            Repeater {
                                model: ui.segmentCount

                                Button {
                                    required property int index
                                    width: 44
                                    height: 32
                                    text: String(index + 1)
                                    checkable: true
                                    checked: ui.selectedSegment === index + 1
                                    highlighted: checked
                                    focusPolicy: Qt.NoFocus
                                    enabled: !ui.graphicsBusy
                                    ButtonGroup.group: steadySegmentButtonGroup
                                    onClicked: Julia.selectSegment(index + 1)
                                }
                            }
                        }
                    }

                    Item {
                        Layout.fillWidth: true
                    }

                    ToolButton {
                        text: "Close"
                        palette.buttonText: "white"
                        onClicked: controlDrawer.close()
                    }
                }
            }

            ScrollView {
                id: steadyValuesScroll
                Layout.fillWidth: true
                Layout.fillHeight: true
                contentHeight: availableHeight
                ScrollBar.vertical.policy: ScrollBar.AlwaysOff
                ScrollBar.horizontal.policy: ScrollBar.AsNeeded
                clip: true

                Row {
                    spacing: 10

                    Item {
                        width: 6
                        height: 1
                    }

                    Repeater {
                        model: window.variables

                        Rectangle {
                            required property int index
                            required property var modelData
                            width: 260
                            height: 54
                            y: Math.max(0, (steadyValuesScroll.availableHeight - height) / 2)
                            radius: 5
                            color: "#e7ebf0"
                            border.color: "#c3cad4"

                            RowLayout {
                                anchors.fill: parent
                                anchors.margins: 6
                                spacing: 7

                                Label {
                                    text: modelData
                                    font.bold: true
                                    Layout.preferredWidth: 42
                                    elide: Text.ElideRight
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
                        width: 6
                        height: 1
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
        anchors.bottom: parent.bottom
        anchors.bottomMargin: Math.max(
            bottomDrawer.visible ? bottomDrawer.height : 0,
            controlDrawer.visible ? controlDrawer.height : 0
        )
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
        sequence: "Space"
        context: Qt.WindowShortcut
        enabled: window.active && !window.textEditorFocused && !ui.graphicsBusy
        autoRepeat: false
        onActivated: Julia.toggleRunning()
    }

    Shortcut {
        sequence: "R"
        context: Qt.WindowShortcut
        enabled: window.active && !window.textEditorFocused && !ui.graphicsBusy
        autoRepeat: false
        onActivated: Julia.resetSimulation()
    }

    Shortcut {
        sequence: "1"
        context: Qt.WindowShortcut
        enabled: window.active && !window.textEditorFocused && !ui.graphicsBusy
        autoRepeat: false
        onActivated: Julia.setDtExponent(-3)
    }

    Shortcut {
        sequence: "2"
        context: Qt.WindowShortcut
        enabled: window.active && !window.textEditorFocused && !ui.graphicsBusy
        autoRepeat: false
        onActivated: Julia.setDtExponent(-1)
    }

    Shortcut {
        sequence: "3"
        context: Qt.WindowShortcut
        enabled: window.active && !window.textEditorFocused && !ui.graphicsBusy
        autoRepeat: false
        onActivated: Julia.setDtExponent(1)
    }

    Shortcut {
        sequence: "4"
        context: Qt.WindowShortcut
        enabled: window.active && !window.textEditorFocused && !ui.graphicsBusy
        autoRepeat: false
        onActivated: Julia.setDtExponent(5)
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
