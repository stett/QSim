#include <SFML/Graphics.hpp>
#include <SFGUI/SFGUI.hpp>
#include "qsim/QSimModel.h"
#include "qsim/QSimView.h"
#include "qsim/QSimController.h"
#include "qsim/WaveFunctionPresets.h"

int main() {

    // Create a window
    sf::RenderWindow window(sf::VideoMode(800, 600), "QSim");

    // Initialize SFGUI
    sfg::SFGUI gui;

    // Initialize QSim stuff
    qsim::QSimModel model(qsim::gaussian_0, qsim::square_well_0);
    qsim::QSimController controller(&model);
    qsim::QSimView view(&model, &controller);

    // While the QSim controller says to keep running, keep runnin'
    while (controller.get_status() == qsim::QSimController::Status::Run) {

        // check all the window's events that were triggered since the last iteration of the loop
        sf::Event event;
        while (window.pollEvent(event)) {
            view.handle_event(event);
        }

        // Update the simulation controller and view
        controller.evolve();
        view.update();

        // Clear the window
        window.clear(sf::Color::Black);

        // Render the view
        view.render();

        // Render the GUI
        gui.Display(window);

        // Display the window
        window.display();
    };

    return 0;
}
