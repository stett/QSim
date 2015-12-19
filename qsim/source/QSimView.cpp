#include <SFML/Graphics.hpp>
#include "qsim/QSimView.h"
#include "qsim/QSimModel.h"

qsim::QSimView::QSimView(QSimModel *model, QSimController *controller) : model(model), controller(controller) {
    window.create(sf::VideoMode(800, 600), "QSim");
}

qsim::QSimView::~QSimView() {
    window.close();
}

void qsim::QSimView::update() {

    // check all the window's events that were triggered since the last iteration of the loop
    sf::Event event;
    while (window.pollEvent(event)) {

        // "close requested" event: we close the window
        if (event.type == sf::Event::Closed)
            controller->quit();
    }
}

void qsim::QSimView::render() {
    window.clear(sf::Color::Black);
    window.display();
}