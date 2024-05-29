"""
 Title:         Model Template
 Description:   Contains the basic structure for a model class
 Author:        Janzen Choi

"""

# Libraries
import importlib, inspect, os, pathlib, sys
import threading

# The Model Template Class
class __Model__:

    def __init__(self, name:str):
        """
        Class for defining a model
        """
        self.name = name
        self.outputs = None

    def get_name(self) -> str:
        """
        Gets the name of the mmodel
        """
        return self.name

    def get_param_names(self) -> list:
        """
        Returns a list of the parameter names for a model
        """
        model_params = inspect.signature(self.run_model).parameters
        param_names = list(model_params.keys())
        return param_names

    def run(self, param_dict:dict, max_time:float=1e5) -> str:
        """
        Calls the run_model function using a thread with timeout;
        also does extra processingg
        
        Parameters:
        * `param_dict`: The parameters for the model as a dictionary
        * `max_time`:   The maximum time to run the model in seconds

        Returns the running status
        """

        # Creates a thread to run the model
        thread = threading.Thread(target=self.__run__, kwargs=param_dict)
        thread.start()
        thread.join(timeout=max_time)

        # Returns the status
        if thread.is_alive():
            return "timeout"
        if self.output == None:
            return "failed"
        return "success"

    def __run__(self, **params) -> tuple:
        """
        Calls the implemented run_model function with a try and stores the outputs
        """
        try:
            output = self.run_model(**params)
            self.output = output
        except:
            self.output = None

    def get_output(self) -> tuple:
        """
        Gets the saved output from the model
        """
        return self.output

    def initialise(self) -> None:
        """
        Runs at the start, once
        """
        pass
        
    def run_model(self, **params) -> tuple:
        """
        Runs the model (must be overridden); returns none if the parameters / model is invalid
        """
        raise NotImplementedError

def create_model(model_name:str, **kwargs) -> __Model__:
    """
    Creates and return a model

    Parameters:
    * `model_name`: The name of the model

    Returns the model
    """

    # Get available models in current folder
    models_dir = pathlib.Path(__file__).parent.resolve()
    files = os.listdir(models_dir)
    files = [file.replace(".py", "") for file in files]
    files = [file for file in files if not file in ["__model__", "__pycache__"]]
    
    # Raise error if model name not in available models
    if not model_name in files:
        raise NotImplementedError(f"The model '{model_name}' has not been implemented")

    # Prepare dynamic import
    module_path = f"{models_dir}/{model_name}.py"
    spec = importlib.util.spec_from_file_location("model_file", module_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    
    # Initialise and return the model
    from model_file import Model
    model = Model(model_name)
    model.initialise(**kwargs)
    return model
