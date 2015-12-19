#ifndef QSIMVIEW_H
#define QSIMVIEW_H
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
        //sf::RenderWindow window;
        //sfg::Window window;
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
        void render();
    };
}

#endif