/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */


#include <cstdlib>     // atexit
// #include <ncurses.h>
// #include <iostream>
// #include <signal.h>     // to handle SIGWINCH signal with handle_winch

#include "curses_funcs.h"


// from http://hughm.cs.ukzn.ac.za/~murrellh/os/notes/ncurses.html
void end_curses(struct CURSES_STRUCT *curses) {
    if (curses->started && !isendwin())
        endwin();
}

void end_fast() {
    endwin();
}

// from http://hughm.cs.ukzn.ac.za/~murrellh/os/notes/ncurses.html
void start_curses(struct CURSES_STRUCT *curses) {
    if (curses->started) {
        refresh();
    } else {
        initscr();
        start_color();
        cbreak();
        noecho();
        std::atexit(end_fast);  // cannot use end_curses with atexit, because atexit
                                // doesn't take arguments. If atexit is omitted, the
                                // terminal will throw a bus error if the program exits
                                // unexpectedly (e.g. due to a bug or error).
        intrflush(stdscr, false);
        keypad(stdscr, true);
        // colors
        init_pair(1, COLOR_BLACK, COLOR_YELLOW); // title
        init_pair(2, COLOR_WHITE, COLOR_BLUE);  // for table
        init_pair(3, COLOR_WHITE, COLOR_BLACK); // for status and message
        // init_pair(4, COLOR_WHITE, COLOR_CYAN);  // for table_fixed 
        curses->started = true;
    }
}

void init_windows(struct CURSES_STRUCT *curses) {
    // initialize windows
    int ncols, nlines;
    getmaxyx(stdscr, nlines, ncols);

    curses->title = newwin(1, ncols, 0, 0);
    curses->message = newwin(1, ncols, 1, 0);
    curses->status = newwin(1, ncols, 2, 0);
    curses->header = newwin(1, ncols, 3, 0);
    // curses->table_fixed = newwin(1, ncols, 2, 0);
    curses->debug = newwin(2, ncols, 4, 0);
    // curses->table = newwin(nlines-4, ncols, 4, 0);  // without debug
    curses->table = newwin(nlines-4-2, ncols, 4+2, 0);  // with debug
    // colors
    wattron(curses->title, COLOR_PAIR(1u));
    wbkgd(curses->title, COLOR_PAIR(1u));
    wattron(curses->table, COLOR_PAIR(2u));
    wbkgd(curses->table, COLOR_PAIR(2u));
    // wattron(curses->table_fixed, COLOR_PAIR(4u));
    // wbkgd(curses->table_fixed, COLOR_PAIR(4u));
    wattron(curses->status, COLOR_PAIR(3u));
    wbkgd(curses->status, COLOR_PAIR(3u));
    wattron(curses->message, COLOR_PAIR(3u));
    wbkgd(curses->message, COLOR_PAIR(3u));
    wattron(curses->header, COLOR_PAIR(2u));
    wbkgd(curses->header, COLOR_PAIR(2u));
    wattron(curses->header, A_BOLD);

    wattron(curses->debug, COLOR_PAIR(1u));
    wbkgd(curses->debug, COLOR_PAIR(1u));

    // scrolling
    scrollok(curses->table, TRUE);
    scrollok(curses->status, TRUE);
    scrollok(curses->message, TRUE);
    scrollok(curses->debug, TRUE);
}
