/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */


#pragma once

#include <ncurses.h>

struct CURSES_STRUCT {
    bool started = false;
    WINDOW *table;
    WINDOW *status;
    WINDOW *title;  // codename + pwd
    WINDOW *header;
    // WINDOW *table_fixed;  // for permanently displaying some rows, e.g. every
                          // order of mag (1, 10, 100, 1000, ...)
    WINDOW *message;
    WINDOW *debug;
};

void end_fast();
void end_curses(struct CURSES_STRUCT *curses);
void start_curses(struct CURSES_STRUCT *curses);
void init_windows(struct CURSES_STRUCT *curses);
