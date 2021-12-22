from pyworkflow.gui import ListTreeProviderString, dialog
from pyworkflow.object import String
from pyworkflow.wizard import Wizard
from prody.protocols import prodyPrefixHelloWorld

class prodyPrefixHelloWorldWizard(Wizard):
    # Dictionary to target protocol parameters
    _targets = [(prodyPrefixHelloWorld, ['message'])]

    def show(self, form, *params):

        # This are the greetings:
        greetings = [String("Hello world"), String("Hola mundo"),
                     String("Bonjour le monde"), String("Hallo Welt"),
                     String("Kon'nichiwa sekai"), String("Nǐ hǎo, shìjiè"),
                     String("Ciao mondo"), String("Hallo Wereld"),
                     String("Privet, mir")]

        # Get a data provider from the greetings to be used in the tree (dialog)
        provider = ListTreeProviderString(greetings)

        # Show the dialog
        dlg = dialog.ListDialog(form.root, "Greetings from the world", provider,
                                "Select one of the greetings)")

        # Set the chosen value back to the form
        form.setVar('message', dlg.values[0].get())
