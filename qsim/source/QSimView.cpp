#include <SFML/Window.hpp>
#include <SFML/Graphics.hpp>
#include <SFGUI/Window.hpp>
#include <SFGUI/Desktop.hpp>
#include "qsim/QSimView.h"
#include "qsim/QSimModel.h"

#define QSIM_COORD_SCREEN_ORIGIN_X(space_min_x, space_range_x, screen_w) (-space_min_x * screen_w / space_range_x)
#define QSIM_COORD_SCREEN_ORIGIN_Y(space_min_y, space_range_y, screen_h) (screen_h * (1.0 + space_min_y / space_range_y))

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

void qsim::QSimView::render(sf::RenderTarget& target, sf::RenderStates states) {

    render_axes(target, states);
}

void qsim::QSimView::render_axes(sf::RenderTarget& target, sf::RenderStates states) {
    double screen_origin_x = QSIM_COORD_SCREEN_ORIGIN_X(model->x_min(), model->x_range(), 800.0);
    double screen_origin_y = QSIM_COORD_SCREEN_ORIGIN_Y(model->y_min(), model->y_range(), 600.0);
    sf::VertexArray vertices(sf::PrimitiveType::Lines, 4);
    vertices[0] = sf::Vertex(sf::Vector2f(screen_origin_x, 0.0f),   sf::Color::White);
    vertices[1] = sf::Vertex(sf::Vector2f(screen_origin_x, 600.0f), sf::Color::White);
    vertices[2] = sf::Vertex(sf::Vector2f(0.0f, screen_origin_y),   sf::Color::White);
    vertices[3] = sf::Vertex(sf::Vector2f(800.0f, screen_origin_y), sf::Color::White);
    target.draw(vertices);
}