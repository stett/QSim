#include "qsim/QSimController.h"
#include "qsim/QSimModel.h"

qsim::QSimController::QSimController(QSimModel *model) : model(model), status(Run) {}

qsim::QSimController::~QSimController() {}

void qsim::QSimController::evolve() {
}

void qsim::QSimController::quit() {
    status = Quit;
}

qsim::QSimController::Status qsim::QSimController::get_status() {
    return status;
}
