#include <string>
#include <SFML/Window.hpp>
#include <SFML/Graphics.hpp>
#include <SFML/System.hpp>
#include <SFGUI/Window.hpp>
#include <SFGUI/Desktop.hpp>
#include <SFGUI/Label.hpp>
#include "qsim/QSimView.h"
#include "qsim/QSimModel.h"
#include "qsim/QSimCoordinates.h"
#include "qsim/complex_colors.h"


qsim::QSimView::QSimView(QSimModel *model, QSimController *controller) : model(model), controller(controller), frames(0.0f), fps(0.0f) {

    // Set up the GUI controls
    auto window = sfg::Window::Create();
    window->SetTitle("Data");

    label_fps = sfg::Label::Create();
    label_normalization = sfg::Label::Create();
    window->Add(label_fps);
    window->Add(label_normalization);

    desktop.Add(window);

    // Set up initial settings
    psi_color = sf::Color(100, 160, 230, 30);
    V_color = sf::Color(255, 180, 100, 255);
    axes_color = sf::Color(80, 80, 80);
    text_color = sf::Color(200, 230, 25);
    tick_color = sf::Color(80, 80, 80);
    ticks_per_pixel = 1.0 / 128.0;
    tick_size = 20.0;
    real_thickness = 3.0;
    imag_thickness = 1.0;
    abs2_thickness = 6.0;
    y_min = -1.0;
    y_max = 2.0;
}

qsim::QSimView::~QSimView() {}

void qsim::QSimView::handle_event(const sf::Event &event) {

    // Let SFGUI handle the event too
    desktop.HandleEvent(event);

    // "close requested" event: we close the window
    if (event.type == sf::Event::Closed)
        controller->quit();

    if (event.type == sf::Event::MouseWheelScrolled) {
        float diff = 0.05f * event.mouseWheelScroll.delta;
        //model->set_x_max(model->x_max() - diff);
        y_max -= diff;
    }
}

void qsim::QSimView::update() {


    // If the time on the clock has reached a second, record the fps
    frames += 1.0f;
    if (clock.getElapsedTime().asSeconds() >= 1.0f) {
        fps = frames / clock.restart().asSeconds();
        frames = 0.0f;
    }

    {   // Update GUI elements
        char str[128];

        sprintf(str, "FPS: %f", fps);
        label_fps->SetText(str);

        sprintf(str, "Norm: %f", model->get_psi_norm());
        label_normalization->SetText(str);
    }

    // Update SFGUI
    desktop.Update(1.0f);
}

void qsim::QSimView::render(sf::RenderTarget& target, const sf::RenderStates &states) {

    // Clear the stage
    target.clear(back_color);

    // Draw the real and imaginary parts of the wavefunction
    render_function_part(target, states, model->get_psi(), psi_color, true);
    render_function_part(target, states, model->get_psi(), psi_color, false);
    render_function_abs2(target, states, model->get_psi(), model->get_psi_abs2());

    // Draw the potential well
    render_function_part(target, states, model->get_V(), V_color, true);
    //render_function_part(target, states, model->get_V(), V_color, false);

    // Draw axes
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
    line.setFillColor(color);
    target.draw(line, states);
}

void qsim::QSimView::render_text(sf::RenderTarget& target, const sf::RenderStates &states, double x, double y, char *str, int size) {
    sf::Text text(str, font, size);
    text.setPosition(sf::Vector2f(x, y));
    text.setColor(text_color);
    target.draw(text, states);
}

void qsim::QSimView::render_axes(sf::RenderTarget& target, const sf::RenderStates &states) {
    double screen_origin_x = QSIM_COORD_SCREEN_ORIGIN_X(model->x_min(), model->x_range(), target.getSize().x);
    double screen_origin_y = QSIM_COORD_SCREEN_ORIGIN_Y(y_min, y_range(), target.getSize().y);
    render_line(target, states, sf::Vector2f(0.0f, screen_origin_y), sf::Vector2f(target.getSize().x, screen_origin_y), 2.0f, axes_color);
    render_line(target, states, sf::Vector2f(screen_origin_x, 0.0f), sf::Vector2f(screen_origin_x, target.getSize().y), 2.0f, axes_color);
}

void qsim::QSimView::render_ticks(sf::RenderTarget& target, const sf::RenderStates &states) {

    // Temp variables
    double origin_x = QSIM_COORD_SCREEN_ORIGIN_X(model->x_min(), model->x_range(), target.getSize().x);
    double origin_y = QSIM_COORD_SCREEN_ORIGIN_Y(y_min, y_range(), target.getSize().y);
    double space_x, space_y;
    double screen_x, screen_y;
    int count_x, count_y;
    char   str[128];

    // Loop through the space-coordinates for the x-axis tick marks, starting at 0 and looping over
    // NOTE: We skip the first one because the text at the origin overlaps with the x-axis
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

        sprintf(str, "%f [m]", space_x);
        render_text(target, states, screen_x + 2.0, screen_y - tick_size, str);

        // Go to the next x-value, looping around if necesarry
        space_x += tick_spacing_x(target.getSize().x);
        if (space_x > model->x_range()) {
            space_x -= model->x_range();
        }
    }

    // Loop through the space-coordinates for the y-axis tick marks, starting at 0 and looping over
    // NOTE: We skip the first one because the text at the origin overlaps with the x-axis
    space_y = tick_spacing_y(target.getSize().y);
    screen_x = origin_x;
    count_y = tick_count_y(target.getSize().y);
    for (int n = 1; n < count_y; ++n) {

        // Get the screen coordinate
        screen_y = QSIM_COORD_SPACE_TO_SCREEN_Y(space_y, origin_y, target.getSize().y, y_range());

        // Draw this tick mark
        render_line(target, states,
            sf::Vector2f(screen_x - tick_size, screen_y),
            sf::Vector2f(screen_x + tick_size, screen_y), 2.0f, tick_color);

        sprintf(str, "%f [MeV]", space_y);
        render_text(target, states, screen_x + tick_size * 0.5f, screen_y + 5, str);

        // Go to the next x-value, looping around if necesarry
        space_y += tick_spacing_y(target.getSize().y);
        if (space_y > y_range()) {
            space_y -= y_range();
        }
    }
}

void qsim::QSimView::render_function_part(sf::RenderTarget& target, const sf::RenderStates &states, const double *f, sf::Color color, bool real) {

    // Temp variables
    static sf::VertexArray vertices(sf::LinesStrip, N);
    //static sf::VertexArray vertices(sf::LinesStrip, 2*N);
    double space_x, space_y;
    double screen_x, screen_y;
    double origin_y = QSIM_COORD_SCREEN_ORIGIN_Y(y_min, y_range(), target.getSize().y);
    double thickness = real ? real_thickness : imag_thickness;

    // Loop through the points
    for (int n = 0; n < N; ++n) {

        // Get this point's x-value in space and screen coordinates
        //space_x = QSIM_COORD_INDEX_TO_SPACE_X(n, model->x_min(), model->x_range());
        screen_x = QSIM_COORD_INDEX_TO_SCREEN_X(n, target.getSize().x);

        // Get this point's y-value in space and screen coordinates
        space_y = real ? GSL_COMPLEX_PACKED_REAL(f, 1, n) : GSL_COMPLEX_PACKED_IMAG(f, 1, n);
        screen_y = QSIM_COORD_SPACE_TO_SCREEN_Y(space_y, origin_y, target.getSize().y, y_range());

        // Set this vertex's properties
        vertices[n].position.x = screen_x;
        vertices[n].position.y = screen_y;
        vertices[n].color = color;
        /*
        int i = 2 * n;
        vertices[i].position.x = screen_x;
        vertices[i].position.y = screen_y;
        vertices[i].color = color;

        vertices[i + 1].position.x = screen_x;
        vertices[i + 1].position.y = screen_y + thickness;
        vertices[i + 1].color = color;
        */
    }

    // Draw the vertices
    target.draw(vertices, states);
}

void qsim::QSimView::render_function_abs2(sf::RenderTarget& target, const sf::RenderStates &states, const double *f, const double *f_abs2) {

    // Temp variables
    //static sf::VertexArray vertices(sf::LinesStrip, N);
    static sf::VertexArray vertices(sf::TrianglesStrip, 2 * N);
    double space_x, space_y;
    double screen_x, screen_y;
    double origin_y = QSIM_COORD_SCREEN_ORIGIN_Y(y_min, y_range(), target.getSize().y);

    // Loop through the points
    for (int n = 0; n < N; ++n) {

        // Get this point's x-value in space and screen coordinates
        screen_x = QSIM_COORD_INDEX_TO_SCREEN_X(n, target.getSize().x);

        // Get this point's y-value in space and screen coordinates
        space_y = f_abs2[n];
        screen_y = QSIM_COORD_SPACE_TO_SCREEN_Y(space_y, origin_y, target.getSize().y, y_range());

        // Get the color representation of this complex value
        sf::Color color = color_complex(GSL_COMPLEX_PACKED_IMAG(f, 1, n));

        // Set this vertex's properties
        int i = 2 * n;
        vertices[i].position.x = screen_x;
        vertices[i].position.y = screen_y;
        vertices[i].color = color;

        vertices[i + 1].position.x = screen_x;
        vertices[i + 1].position.y = screen_y + abs2_thickness;
        vertices[i + 1].color = color;
    }

    // Draw the vertices
    target.draw(vertices, states);
}

void qsim::QSimView::set_font(const std::string &file_path) {
    font.loadFromFile(file_path);
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
    return y_range() / (screen_h * ticks_per_pixel);
}

double qsim::QSimView::y_range() {
    return y_max - y_min;
}
