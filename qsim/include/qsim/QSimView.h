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
        void render(sf::RenderTarget& target, sf::RenderStates states);

        // Internal draw methods
    private:
        void render_axes(sf::RenderTarget& target, sf::RenderStates states);
    };
}

#endif