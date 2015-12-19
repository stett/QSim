#include <SFML/Window.hpp>
#include <SFML/Graphics.hpp>
#include <SFGUI/Window.hpp>
#include <SFGUI/Desktop.hpp>
#include "qsim/QSimView.h"
#include "qsim/QSimModel.h"
#include "qsim/QSimCoordinates.h"


qsim::QSimView::QSimView(QSimModel *model, QSimController *controller) : model(model), controller(controller) {
    auto window = sfg::Window::Create();
    window->SetTitle("Test Window");
    desktop.Add(window);

    // Set up initial settings
    axes_color = sf::Color(100, 100, 100);
    text_color = sf::Color(200, 230, 25);
    tick_color = sf::Color(230, 230, 230);
    ticks_per_pixel = 1.0 / 128.0;
    tick_size = 20.0;

    // Build meshes
    vertices = sf::VertexArray(sf::PrimitiveType::Lines, 4);
    vertices[0] = sf::Vertex(sf::Vector2f( 0.0f, -1.0f), axes_color);
    vertices[1] = sf::Vertex(sf::Vector2f( 0.0f,  1.0f), axes_color);
    vertices[2] = sf::Vertex(sf::Vector2f(-1.0f,  0.0f), axes_color);
    vertices[3] = sf::Vertex(sf::Vector2f( 1.0f,  0.0f), axes_color);
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

void qsim::QSimView::render(sf::RenderTarget& target, const sf::RenderStates &states) {
    target.clear(back_color);
    render_axes(target, states);
    render_ticks(target, states);
}

void qsim::QSimView::render_line(sf::RenderTarget& target, const sf::RenderStates &states, const sf::Vector2f &p0, const sf::Vector2f &p1, float thickness, sf::Color color) {
    float dx = p1.x - p0.x;
    float dy = p1.y - p0.y;
    float dist = sqrt(dx*dx + dy*dy);
    float angle = atan2(dy, dx);
    sf::RectangleShape line(sf::Vector2f(dist, thickness));
    line.setPosition(p0);
    line.rotate(angle * 180.0f / (float)M_PI);
    line.setOutlineColor(color);
    target.draw(line, states);
}

void qsim::QSimView::render_axes(sf::RenderTarget& target, const sf::RenderStates &states) {
    double screen_origin_x = QSIM_COORD_SCREEN_ORIGIN_X(model->x_min(), model->x_range(), target.getSize().x);
    double screen_origin_y = QSIM_COORD_SCREEN_ORIGIN_Y(model->y_min(), model->y_range(), target.getSize().y);
    render_line(target, states, sf::Vector2f(0.0f, screen_origin_y), sf::Vector2f(target.getSize().x, screen_origin_y), 2.0f, axes_color);
    render_line(target, states, sf::Vector2f(screen_origin_x, 0.0f), sf::Vector2f(screen_origin_x, target.getSize().y), 2.0f, axes_color);
}

void qsim::QSimView::render_ticks(sf::RenderTarget& target, const sf::RenderStates &states) {
    //glLineWidth(1.f);
    //tick_color.setGlColor();

    // Temp variables
    double origin_x = QSIM_COORD_SCREEN_ORIGIN_X(model->x_min(), model->x_range(), target.getSize().x);
    double origin_y = QSIM_COORD_SCREEN_ORIGIN_Y(model->y_min(), model->y_range(), target.getSize().y);
    double space_x, space_y;
    double screen_x, screen_y;
    int count_x, count_y;
    char   str[10];

    // Loop through the space-coordinates for the x-axis tick marks, starting at 0 and looping over
    // NOTE: We skip the first one because the text at the origin overlaps with the x-axis
    //    space_x  = 0;
    space_x = tick_spacing_x(target.getSize().x);
    screen_y = origin_y;
    count_x = tick_count_x(target.getSize().x);
    for (int n = 1; n < count_x; ++n) {

        // Get the screen coordinate
        screen_x = QSIM_COORD_SPACE_TO_SCREEN_X(space_x, origin_x, target.getSize().x, model->x_range());

        // Draw this tick mark
        render_line(target, states,
            sf::Vector2f(screen_x, screen_y - tick_size),
            sf::Vector2f(screen_x, screen_y + tick_size), 2.0f, tick_color);

        //sprintf(str, "%f", space_x);
        //draw_string(screen_x + 2., screen_y - .5 * tick_size, str);

        // Go to the next x-value, looping around if necesarry
        space_x += tick_spacing_x(target.getSize().x);
        if (space_x > model->x_range()) {
            space_x -= model->x_range();
        }
    }

    // Loop through the space-coordinates for the y-axis tick marks, starting at 0 and looping over
    // NOTE: We skip the first one because the text at the origin overlaps with the x-axis
    //    space_y  = 0;
    space_y = tick_spacing_y(target.getSize().y);
    screen_x = origin_x;
    count_y = tick_count_y(target.getSize().y);
    for (int n = 1; n < count_y; ++n) {

        // Get the screen coordinate
        screen_y = QSIM_COORD_SPACE_TO_SCREEN_Y(space_y, origin_y, target.getSize().y, model->y_range());

        // Draw this tick mark
        render_line(target, states,
            sf::Vector2f(screen_x - tick_size, screen_y),
            sf::Vector2f(screen_x + tick_size, screen_y), 2.0f, tick_color);

        //sprintf(str, "%f", space_y);
        //draw_string(screen_x + tick_size + 1, screen_y + 5, str);

        // Go to the next x-value, looping around if necesarry
        space_y += tick_spacing_y(target.getSize().y);
        if (space_y > model->y_range()) {
            space_y -= model->y_range();
        }
    }
}

double qsim::QSimView::tick_count_x(double screen_w) {
    return screen_w * ticks_per_pixel;
}

double qsim::QSimView::tick_count_y(double screen_h) {
    return screen_h * ticks_per_pixel;
}

double qsim::QSimView::tick_spacing_x(double screen_w) {
    return model->x_range() / (screen_w * ticks_per_pixel);
}

double qsim::QSimView::tick_spacing_y(double screen_h) {
    return model->y_range() / (screen_h * ticks_per_pixel);
}
