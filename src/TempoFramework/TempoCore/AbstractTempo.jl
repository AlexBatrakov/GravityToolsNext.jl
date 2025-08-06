#--------------------------------------------------------------------------------------------------------------
# Abstract type for different versions of Tempo.
"""
    AbstractTempoVersion

An abstract type that serves as a base for different implementations of Tempo versions.
Implementations must define the `get_tempo_command` function.
"""
abstract type AbstractTempoVersion end

# Function definition for all Tempo types
"""
    get_tempo_command(tempo_version::AbstractTempoVersion)

Obtain the command line call string for the specified Tempo version. This function must be
implemented for each specific subtype of `AbstractTempoVersion`.
"""
function get_tempo_command(::AbstractTempoVersion)
    error("This function must be implemented for each specific Tempo version.")
end

# Constructors and functions for Tempo
"""
    Tempo(custom_data_directory::String = get_tempo_directory())

A structure representing the original Tempo version with an optional custom data directory.
If no directory is provided, the default is retrieved from the `TEMPO` environment variable.

# Constructors
- `Tempo()`: Uses the default directory path from `TEMPO` environment variable.
- `Tempo(custom_dir::String)`: Uses a custom directory specified by `custom_dir`.
"""
struct Tempo <: AbstractTempoVersion
    custom_data_directory::String
    Tempo() = new(get_tempo_directory()) # Default constructor
    Tempo(custom_dir::String) = new(custom_dir) # Constructor with a parameter
end

# Function to get the default directory for Tempo
"""
    get_tempo_directory() -> String

Retrieve the default directory path for Tempo from the `TEMPO` environment variable.
If the variable is not set, it returns a predefined default path.
"""
get_tempo_directory() = get(ENV, "TEMPO", "/default/path/to/tempo")

# Constructors and functions for Tempo2
"""
    Tempo2(custom_data_directory::String = get_tempo2_directory())

A structure representing the Tempo2 version with an optional custom data directory.
If no directory is provided, the default is retrieved from the `TEMPO2` environment variable.

# Constructors
- `Tempo2()`: Uses the default directory path from `TEMPO2` environment variable.
- `Tempo2(custom_dir::String)`: Uses a custom directory specified by `custom_dir`.
"""
struct Tempo2 <: AbstractTempoVersion
    custom_data_directory::String
    Tempo2() = new(get_tempo2_directory()) # Default constructor
    Tempo2(custom_dir::String) = new(custom_dir) # Constructor with a parameter
end

# Function to get the default directory for Tempo2
"""
    get_tempo2_directory() -> String

Retrieve the default directory path for Tempo2 from the `TEMPO2` environment variable.
If the variable is not set, it returns a predefined default path.
"""
get_tempo2_directory() = get(ENV, "TEMPO2", "/default/path/to/tempo2")

# Function definitions for getting the command
"""
    get_tempo_command(tempo_version::Tempo) -> String

Returns the command line call string for the original Tempo version.
"""
get_tempo_command(::Tempo) = "tempo"

"""
    get_tempo_command(tempo_version::Tempo2) -> String

Returns the command line call string for the Tempo2 version.
"""
get_tempo_command(::Tempo2) = "tempo2"

