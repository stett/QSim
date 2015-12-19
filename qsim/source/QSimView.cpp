#include <SFML/Window.hpp>
#include <SFML/Graphics.hpp>
#include <SFGUI/Window.hpp>
#include <SFGUI/Desktop.hpp>
#include "qsim/QSimView.h"
#include "qsim/QSimModel.h"

qsim::QSimView::QSimView(QSimModel *model, QSimController *controller) : model(model), controller(controller) {
    auto window = sfg::Window::Create();
    window->SetTitle("Test Window");
    desktop.Add(window);
}

qsim::QSimView::~QSimView() {}

void qsim::QSimView::handle_event(const sf::Event &event) {

    // Let SFGUI handle the event too
    desktop.HandleEvent(event);

    // "close requested" event: we close the window
    if (event.type == sf::Event::Closed)
        controller->quit();
}

void qsim::QSimView::update() {

    // Update SFGUI
    desktop.Update(1.0f);
}

void qsim::QSimView::render() {
}
