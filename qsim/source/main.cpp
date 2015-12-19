#include <SFGUI/SFGUI.hpp>
#include "qsim/QSimModel.h"
#include "qsim/QSimView.h"
#include "qsim/QSimController.h"
#include "qsim/WaveFunctionPresets.h"

int main() {

    // Initialize SFGUI
    sfg::SFGUI sfgui;

    // Initialize QSim stuff
    qsim::QSimModel model(qsim::gaussian_0, qsim::square_well_0);
    qsim::QSimController controller(&model);
    qsim::QSimView view(&model, &controller);

    // While the QSim controller says to keep running, keep runnin'
    while (controller.get_status() == qsim::QSimController::Status::Run) {
        controller.evolve();
        view.update();
        view.render();
    };

    return 0;
}
