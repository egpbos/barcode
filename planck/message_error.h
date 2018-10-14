/*
 *  Class for error reporting
 *
 *  Copyright (C) 2003, 2004 Max-Planck-Society
 *  Authors: Reinhard Hell, Martin Reinecke
 */

#ifndef PLANCK_MESSAGE_ERROR_H
#define PLANCK_MESSAGE_ERROR_H

#include <exception>
#include <iostream>
#include <string>

#if defined (PLANCK_STACKTRACE)
#include <execinfo.h>
#endif

inline void show_stackframe()
  {
#if defined (PLANCK_STACKTRACE)
  void *trace[16];
  int trace_size = backtrace(trace, 16);
  char **messages = backtrace_symbols(trace, trace_size);
  std::cerr << "[bt] Execution path:" << std::endl;
  for (int i=0; i<trace_size; ++i)
    std::cerr << "[bt] " << messages[i] << std::endl;
#endif
  }


class Message_error
  {
  private:
    std::string msg;

  public:
    Message_error()
      : msg (std::string("Unspecified error"))
      { std::cerr<<msg<<std::endl; show_stackframe(); }

    explicit Message_error(const std::string &message)
      : msg (message) { std::cerr<<msg<<std::endl; show_stackframe(); }

    virtual const char* what() const
      { return msg.c_str(); }

    virtual ~Message_error() {}
  };

#if defined (PLANCK_CHECKS)

#define PLANCK_DIAGNOSIS_BEGIN try {
#define PLANCK_DIAGNOSIS_END \
} \
catch (Message_error &e) \
  { std::cerr << "Planck exception: " << e.what() << std::endl; throw; } \
catch (std::exception &e) \
  { std::cerr << "std::exception: " << e.what() << std::endl; throw; } \
catch (...) \
  { std::cerr << "Unknown exception" << std::endl; throw; }

#else

#define PLANCK_DIAGNOSIS_BEGIN
#define PLANCK_DIAGNOSIS_END

#endif

#endif
