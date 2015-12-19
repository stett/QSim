#ifndef QSIMVIEW_H
#define QSIMVIEW_H
#include <SFML/Graphics.hpp>
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
        sf::RenderWindow window;

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
        void update();
        void render();
    };
}

#endif