#ifndef QSIMVIEW_H
#define QSIMVIEW_H
#include <memory>
#include <string>
#include <SFML/Graphics.hpp>
#include <SFML/System.hpp>
#include <SFGUI/Window.hpp>
#include <SFGUI/Desktop.hpp>
#include <SFGUI/Label.hpp>
#include <SFGUI/Box.hpp>
#include "qsim/QSimModel.h"
#include "qsim/QSimController.h"

namespace qsim {
    class QSimView {

        // MVC data
    private:
        QSimModel *model;
        QSimController *controller;

        // View data
    private:
        sfg::Desktop desktop;
        sf::Font font;
        sf::Clock clock;
        float frames;
        float fps;

        // GUI data
    private:
        std::shared_ptr<sfg::Box> data_box;
        std::shared_ptr<sfg::Label> fps_label;
        std::shared_ptr<sfg::Label> norm_label;

        // Settings
    private:
        sf::Color back_color, axes_color, text_color, tick_color;
        sf::Color psi_color, V_color;
        double real_thickness, imag_thickness, abs2_thickness;
        double ticks_per_pixel;
        double tick_size;
        double y_min, y_max;

        // 'Tors
    public:
        QSimView(QSimModel *model, QSimController *controller);
        ~QSimView();

        // Non-copyable
    public:
        QSimView(const QSimView &view) = delete;
        QSimView &operator=(const QSimView &view) = delete;

        // Public methods
    public:
        void handle_event(const sf::Event &event);
        void update();
        void render(sf::RenderTarget& target, const sf::RenderStates &states);

        // Internal draw methods
    private:
        void render_line(sf::RenderTarget& target, const sf::RenderStates &states, const sf::Vector2f &p0, const sf::Vector2f &p1, float thickness = 1.0f, sf::Color color = sf::Color::White);
        void render_text(sf::RenderTarget& target, const sf::RenderStates &states, double x, double y, char *str, int size = 12);
        void render_axes(sf::RenderTarget& target, const sf::RenderStates &states);
        void render_ticks(sf::RenderTarget& target, const sf::RenderStates &states);
        void render_function_part(sf::RenderTarget& target, const sf::RenderStates &states, const double *f, sf::Color colof = sf::Color::White, bool real = true);
        void render_function_abs2(sf::RenderTarget& target, const sf::RenderStates &states, const double *f, const double *f_abs2);

        // Setters
    public:
        void set_font(const std::string &file_path);

        // Getters
    public:
        double tick_count_x(double screen_w);
        double tick_count_y(double screen_h);
        double tick_spacing_x(double screen_w);
        double tick_spacing_y(double screen_h);
        double y_range();
    };
}

#endif