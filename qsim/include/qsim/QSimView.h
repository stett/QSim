#ifndef QSIMVIEW_H
#define QSIMVIEW_H
#include <SFML/Graphics.hpp>
#include <SFGUI/Window.hpp>
#include <SFGUI/Desktop.hpp>
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

        // Meshes
    private:
        sf::VertexArray vertices;

        // Settings
    private:
        sf::Color back_color, axes_color, text_color, tick_color;
        double ticks_per_pixel;
        double tick_size;

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
        void render_axes(sf::RenderTarget& target, const sf::RenderStates &states);
        void render_ticks(sf::RenderTarget& target, const sf::RenderStates &states);

        // Getters
    public:
        double tick_count_x(double screen_w);
        double tick_count_y(double screen_h);
        double tick_spacing_x(double screen_w);
        double tick_spacing_y(double screen_h);
    };
}

#endif