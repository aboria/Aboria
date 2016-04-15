from jinja2 import Environment, FileSystemLoader
import os

# Capture our current directory
THIS_DIR = os.path.dirname(os.path.abspath(__file__))

class Variable:
    def __init__(self, type_name, representation, description):
        self.type_name = type_name
        self.representation = representation
        self.description = description

class Particles:
    def __init__(self,list_of_variables):
        self.name = "particles"
        self.particle_name = "particle"
        list_of_names = [variable.type_name for variable in list_of_variables]
        list_of_names_referenced = [variable.type_name+'&' for variable in list_of_variables]

        # Create the jinja2 environment.
        # Notice the use of trim_blocks, which greatly helps control whitespace.
        j2_env = Environment(loader=FileSystemLoader(THIS_DIR),
                        block_start_string='/*{%',
                        block_end_string='%}*/',
                        variable_start_string='/*{{',
                        variable_end_string='}}*/',
                        trim_blocks=True
                        )
        output_from_parsed_template = j2_env.get_template('python_template.cpp').render(
            variable_type_string = 'std::tuple<' + ', '.join(list_of_names) + '>',
            value_type_string = 'std::tuple<' + ', '.join(list_of_names_referenced) + '>',
            particles = self,
            variables = list_of_variables
        )

        with open('%s.cpp'%self.name, "wb") as f:
            f.write(output_from_parsed_template)


v = Variable('velocity','Vect3d','this is the velocity')
p = Variable('pressure','double','this is the pressure')

particles = Particles([v,p])
