#ifndef SEMAPP_H
#define SEMAPP_H

#include "MooseApp.h"

class SemApp;

template<>
InputParameters validParams<SemApp>();

class SemApp : public MooseApp
{
public:
  SemApp(const std::string & name, InputParameters parameters);
  virtual ~SemApp();

  static void registerApps();
  static void registerObjects(Factory & factory);
  static void associateSyntax(Syntax & syntax, ActionFactory & action_factory);
};

#endif /* SEMAPP_H */
