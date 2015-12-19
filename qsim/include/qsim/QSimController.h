#ifndef QSIMCONTROLLER_H
#define QSIMCONTROLLER_H
#include "qsim/QSimModel.h"

namespace qsim {
    class QSimController {

    public:
        enum Status { Run, Quit };

        // Private members
    private:
        QSimModel *model;
        Status status;

        // 'Tors
    public:
        QSimController(QSimModel *model);
        ~QSimController();

        // Non-copyable
    public:
        QSimController(const QSimController &controller) = delete;
        QSimController &operator=(const QSimController &controller) = delete;

        // Public methods
    public:
        void evolve();
        void quit();
        Status get_status();
    };
}

#endif