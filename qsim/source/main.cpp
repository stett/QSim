
#include <stdio.h>
#include <gtest/gtest.h>
#include "test_QSimMath.h"

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    int err = RUN_ALL_TESTS();
    char str = ' ';
    std::cin >> str;
    return err;
}



/*
#include <string>
#include <SFML/Graphics.hpp>
#include <SFGUI/SFGUI.hpp>
#include "qsim/QSimModel.h"
#include "qsim/QSimView.h"
#include "qsim/QSimController.h"
#include "qsim/WaveFunctionPresets.h"


int main(int argc, char **argv) {

    // Load the font for the view
    // Get the last position of '/'
    std::string aux(argv[0]);

    // get '/' or '\\' depending on unix/mac or windows.
#if defined(_WIN32) || defined(WIN32)
    int pos = aux.rfind('\\');
#else
    int pos = aux.rfind('/');
#endif

    // Get the path and the name
    std::string path = aux.substr(0, pos + 1);
    std::string name = aux.substr(pos + 1);

    // Create a window
    sf::ContextSettings settings;
    settings.antialiasingLevel = 4;
    settings.majorVersion = 3;
    settings.minorVersion = 1;
    sf::RenderWindow window(sf::VideoMode(800, 600), "QSim", sf::Style::Default, settings);

    // Initialize SFGUI
    sfg::SFGUI gui;

    //
    qsim::QSimModel model;
    qsim::QSimController controller(&model);
    qsim::QSimView view(&model, &controller);

    // Initialize QSim stuff
    double x_min        = 0.0;
    double x_max        = 100.0;
    double x0           = (x_min + x_max) * 0.3;
    double k0           = 10.0;
    double alpha        = 1.0;
    double V_max        = 8.0;
    double thickness    = 0.48 * (x_max - x_min);
    qsim::Gaussian psi_0(x0, k0, alpha);
    qsim::SquareBarrier V_0(x_min + thickness, x_max - thickness, V_max);
    model.set_functions(psi_0, V_0);

    // Load the font for the view
    view.set_font(path + std::string("ubuntu_mono_regular.ttf"));

    // While the QSim controller says to keep running, keep runnin'
    while (controller.get_status() == qsim::QSimController::Status::Run) {

        // Set up render states
        sf::RenderStates states;

        // check all the window's events that were triggered since the last iteration of the loop
        sf::Event event;
        while (window.pollEvent(event)) {

            // Let the view handle the event
            view.handle_event(event);

            // Resize the viewmatrix
            if (event.type == sf::Event::Resized) {
                sf::Vector2f size(window.getSize().x, window.getSize().y);
                sf::Vector2f center(size.x * 0.5f, size.y * 0.5f);
                window.setView(sf::View(center, size));
            }
        }

        // Update the simulation controller and view
        controller.update();
        view.update();

        // Clear the window
        window.clear(sf::Color::Black);

        // Render the view
        view.render(window, states);

        // Render the GUI
        gui.Display(window);

        // Display the window
        window.display();
    };

    return 0;
}
*/
