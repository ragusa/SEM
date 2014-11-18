#include "SemApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"

template<>
InputParameters validParams<SemApp>()
{
  InputParameters params = validParams<MooseApp>();
  return params;
}

SemApp::SemApp(const std::string & name, InputParameters parameters) :
    MooseApp(name, parameters)
{
  srand(processor_id());

  Moose::registerObjects(_factory);
  ModulesApp::registerObjects(_factory);
  SemApp::registerObjects(_factory);

  Moose::associateSyntax(_syntax, _action_factory);
  ModulesApp::associateSyntax(_syntax, _action_factory);
  SemApp::associateSyntax(_syntax, _action_factory);
}

SemApp::~SemApp()
{
}

void
SemApp::registerApps()
{
  registerApp(SemApp);
}

void
SemApp::registerObjects(Factory & factory)
{
}

void
SemApp::associateSyntax(Syntax & syntax, ActionFactory & action_factory)
{
}
