# Overview
This IFC to gbXML converter is a result of research performed during my graduation project at the Eindhoven University of Technology. The converter aims to create an improved interoperability between BIM and whole-building energy analysis by translating the open exchange format IFC to a validated gbXML file format. Relationships between gbXML elements and corresponding IFC entities are created and linkages are established by iterating over multiple IFC entities.

* Author: Maarten Visschers
* Refactoring and IFC4 support: Bruno Postle

## Using the converter
The current version of the converter is Python based without a GUI, with the script you are able to run the IFC-to-gbXML conversion. Supported are:
* The IFC2X3 and IFC4 schema;
* The gbXML 6.01 schema (http://www.gbxml.org/Schema_Current_GreenBuildingXML_gbXML).

### Command-line usage:
    ifcgbxml.py input.ifc output.xml

### Python usage:
    import ifcopenshell
    from ifcgbxml import create_gbxml
    ifc_file = ifcopenshell.open("/path/input.ifc")
    root = create_gbxml(ifc_file)
    root.writexml(open("/path/output.xml", "w"), indent="  ", addindent="  ", newl="\n")

## What to do
The converter expects that the IFC file has sufficient data, ie:

* Each *Site* should have location attributes, ie. *RefLatitude*, *RefLongitude*, *RefElevation* and *SiteAddress*. *RefLatitude* **must** be correct.
* Each *Building* should have a *BuildingAddress*.
* *Spatial Elements* **must** be structured correctly, ie. a *Space* **must** be contained within a *Building Storey*, which **must** be contained within a *Building*, which **must** be contained within a *Site*.
* Each *Space* **should** have *Qto_SpaceBaseQuantities/NetFloorArea* and *Qto_SpaceBaseQuantities/NetVolume* quantities.
* Each *Space* **must** be bounded by *2nd Level* *Space Boundary* relationships.
* Each *Space Boundary* **must** have a *Related Building Element*, planar *Connection Geometry*, and an *Internal Or External Boundary* attribute.
* Each *Building Element* related to a *Space Boundary* must be a *Window* or a *Door*, or have a *Material Layer Set* construction (eg. a *Wall*, *Roof*, or a *Slab*).
* Each *Window* and *Door* (or its Type) **must** have a U-value property in *Pset_WindowCommon/ThermalTransmittance* or *Pset_DoorCommon/ThermalTransmittance*.
* Each *Window* and *Door* should have a transmittance property in *Pset_WindowCommon/GlazingAreaFraction* or *Pset_DoorCommon/GlazingAreaFraction*.
* Each *Building Element* with a *Material Layer Set* construction (eg. a *Wall, *Roof*, or a *Slab) **must** have a U-value property (eg. *Pset_WallCommon/ThermalTransmittance*) **or** be constructed from layers with full thermal properties.
* Each *Material Layer* in a *Material Layer Set* construction **must** have a *Thickness* attribute, and each *Material* in a *Material Layer* **must** have an R-value property in *Pset_MaterialEnergy/ThermalConductivityTemperatureDerivative*.
* For metric projects, units for U-Values **must** be `W/m²·K`, and R-values **must** be `m²·K/W`. For imperial projects (ie. where the LENGTHUNIT is feet or inches), U-value units **must** be `Btu/h·ft²·F`, and R-values **must** be `h·ft²·F/Btu`.

## Capabilities
* Process basic IFC building elements (e.g. walls and floors);
* Convert implicit to explicit geometric information (to comply with gbXML);
* Process IFC window and opening elements;
* Write valid information according to the official gbXML schema;
* Create gbXML files which are suitable for whole-building energy simulation in DesignBuilder;
* Include gbXML construction types with each their corresponding layer and materials;
* Include thermal properties (e.g. IFC analytical property sets);
* gbXML output is in SI units, project units are converted to SI units.

## More information
For more information, please contact maartenvisschers@hotmail.com, or file an issue in the repository.
