#ifndef QSIMCOORDINATES_H
#define QSIMCOORDINATES_H
#include "qsim/QSimConstants.h"

#define QSIM_COORD_INDEX_TO_SPACE_X(n, space_min_x, space_range_x) ((double)n * space_range_x / (double)N + space_min_x)
#define QSIM_COORD_INDEX_TO_SCREEN_X(n, screen_w) ((double)n * screen_w / (double)(N - 1.0))
#define QSIM_COORD_SCREEN_TO_SPACE_X(screen_x, screen_origin_x, screen_w, space_range_x) ((screen_x - screen_origin_x) * space_range_x / screen_w)
#define QSIM_COORD_SCREEN_TO_SPACE_Y(screen_y, screen_origin_y, screen_h, space_range_y) ((screen_origin_y - screen_y) * space_range_y / screen_h)
#define QSIM_COORD_SPACE_TO_SCREEN_X(space_x, screen_origin_x, screen_w, space_range_x) ( space_x * screen_w / space_range_x + screen_origin_x)
#define QSIM_COORD_SPACE_TO_SCREEN_Y(space_y, screen_origin_y, screen_h, space_range_y) (-space_y * screen_h / space_range_y + screen_origin_y)
#define QSIM_COORD_SCREEN_ORIGIN_X(space_min_x, space_range_x, screen_w) (-space_min_x * screen_w / space_range_x)
#define QSIM_COORD_SCREEN_ORIGIN_Y(space_min_y, space_range_y, screen_h) (screen_h * (1.0 + space_min_y / space_range_y))

#endif