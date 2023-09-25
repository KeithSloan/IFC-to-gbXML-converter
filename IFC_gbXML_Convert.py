#!/usr/bin/python3
# Import necessary python libraries e.g. IfcOpenShell and MiniDom
import ifcopenshell.util.placement
import numpy as np
import datetime
import time
import sys
from xml.dom import minidom


# from https://github.com/wassimj/TopologicSverchok/blob/main/nodes/Topologic/ifc_topologic.py
def get_boundary_vertices(ifc_rel_space_boundary):
    ifc_curve = ifc_rel_space_boundary.ConnectionGeometry.SurfaceOnRelatingElement
    if ifc_curve.is_a("IfcFaceSurface"):
        ifc_points = ifc_curve.Bounds[0].Bound.Polygon
        ifc_plane = ifc_curve.FaceSurface
    elif ifc_curve.is_a("IfcCurveBoundedPlane"):
        ifc_points = ifc_curve.OuterBoundary.Points
        ifc_plane = ifc_curve.BasisSurface
    plane_matrix = ifcopenshell.util.placement.get_axis2placement(ifc_plane.Position)
    plane_coords = [v.Coordinates for v in ifc_points]
    plane_vertices = [
        np.array(v + (0, 1)) if len(v) == 2 else np.array(v + (1,))
        for v in plane_coords
    ]
    vertices = [(plane_matrix @ v)[0:3] for v in plane_vertices]
    if (
        np.dot(
            np.cross(vertices[1] - vertices[0], vertices[-1] - vertices[0]),
            (plane_matrix @ np.array([0, 0, 1, 0]))[0:3],
        )
        < -1e-6
    ):
        vertices.reverse()
    return vertices


def get_poly_loop(root, vertices, ifc_relating_space):
    poly_loop = root.createElement("PolyLoop")
    new_z = ifc_relating_space.ObjectPlacement.PlacementRelTo.RelativePlacement.Location.Coordinates[
        2
    ]
    for v in vertices:
        x, y, z = v
        z = z + new_z
        cartesian_point = root.createElement("CartesianPoint")
        for c in x, y, z:
            if "e-" in str(c):
                c = 0.0
            coordinate = root.createElement("Coordinate")
            coordinate.appendChild(root.createTextNode(str(c)))
            cartesian_point.appendChild(coordinate)
        poly_loop.appendChild(cartesian_point)
    return poly_loop


# Align the gbXML input according to the predefined official gbXML schema
def remove_unnecessary_characters(element):
    char_to_replace = {"$": "", ":": "", " ": "", "(": "", ")": ""}
    for key, value in char_to_replace.items():
        element = element.replace(key, value)
    return element


def fix_xml_cmps(element):
    return "campus" + remove_unnecessary_characters(element)


def fix_xml_bldng(element):
    return "building" + remove_unnecessary_characters(element)


def fix_xml_stry(element):
    return "storey" + remove_unnecessary_characters(element)


def fix_xml_spc(element):
    return "space" + remove_unnecessary_characters(element)


def fix_xml_id(element):
    return "id" + remove_unnecessary_characters(element)


def fix_xml_name(element):
    return "object" + remove_unnecessary_characters(element)


def fix_xml_cons(element):
    return "construct" + remove_unnecessary_characters(element)


def fix_xml_layer(element):
    return "lyr" + remove_unnecessary_characters(element)


if not len(sys.argv) == 3:
    sys.exit("Usage: " + sys.argv[0] + " input.ifc output.xml")
ifc_file = ifcopenshell.open(sys.argv[1])

# Create the XML root by making use of MiniDom
root = minidom.Document()

# Create the 'gbXML' element and append it to the Root of the document
gbxml = root.createElement("gbXML")
root.appendChild(gbxml)

# Create attributes for the 'gbXML' element
gbxml.setAttribute("xmlns", "http://www.gbxml.org/schema")
gbxml.setAttribute("temperatureUnit", "C")
gbxml.setAttribute("lengthUnit", "Meters")
gbxml.setAttribute("areaUnit", "SquareMeters")
gbxml.setAttribute("volumeUnit", "CubicMeters")
gbxml.setAttribute("useSIUnitsForResults", "true")
gbxml.setAttribute("version", "0.37")

# Create a dictionary to store all gbXML element Ids
dict_id = {}

# Specify the 'Campus' element of the gbXML schema; making use of IFC entity 'IfcSite'
# This element is added as child to the earlier created 'gbXML' element
# Site must have the Longitude attribute

ifc_sites = ifc_file.by_type("IfcSite")
ifc_sites = [ifc_site for ifc_site in ifc_sites if ifc_site.RefLongitude != None]

for ifc_site in ifc_sites:
    campus = root.createElement("Campus")
    campus.setAttribute("id", fix_xml_cmps(ifc_site.GlobalId))
    gbxml.appendChild(campus)

    dict_id[fix_xml_cmps(ifc_site.GlobalId)] = campus

    # Specify the 'Location' element of the gbXML schema; making use of IFC entities 'IfcSite' and 'IfcPostalAddress'
    # This new element is added as child to the earlier created 'Campus' element
    location = root.createElement("Location")
    campus.appendChild(location)

    longitude = root.createElement("Longitude")
    longitudeValue = str(ifc_site.RefLongitude[0])
    longitude.appendChild(root.createTextNode(longitudeValue))
    location.appendChild(longitude)

    latitude = root.createElement("Latitude")
    latitudeValue = str(ifc_site.RefLatitude[0])
    latitude.appendChild(root.createTextNode(latitudeValue))
    location.appendChild(latitude)

    elevation = root.createElement("Elevation")
    elevation.appendChild(root.createTextNode(str(ifc_site.RefElevation)))
    location.appendChild(elevation)

    # FIXME assigns all Postal Addresses to all Locations
    for ifc_postal_address in ifc_file.by_type("IfcPostalAddress"):
        zipcode = root.createElement("ZipcodeOrPostalCode")
        zipcode.appendChild(root.createTextNode(ifc_postal_address.PostalCode))
        location.appendChild(zipcode)

        name = root.createElement("Name")
        name.appendChild(
            root.createTextNode(
                ifc_postal_address.Region + ", " + ifc_postal_address.Country
            )
        )
        location.appendChild(name)

    # Specify the 'Building' element of the gbXML schema; making use of IFC entity 'IfcBuilding'
    # This new element is added as child to the earlier created 'Campus' element
    # FIXME assigns all Buildings to all Campuses
    for ifc_building in ifc_file.by_type("IfcBuilding"):
        building = root.createElement("Building")
        building.setAttribute("id", fix_xml_bldng(ifc_building.GlobalId))
        building.setAttribute("buildingType", "Unknown")
        campus.appendChild(building)

        dict_id[fix_xml_bldng(ifc_building.GlobalId)] = building

        # FIXME assigns all Postal Addresses to all buildings
        for ifc_postal_address in ifc_file.by_type("IfcPostalAddress"):
            street_address = root.createElement("StreetAddress")
            street_address.appendChild(
                root.createTextNode(
                    ifc_postal_address.Region + ", " + ifc_postal_address.Country
                )
            )
            building.appendChild(street_address)

        # Specify the 'BuildingStorey' element of the gbXML schema; making use of IFC entity 'IfcBuildingStorey'
        # This new element is added as child to the earlier created 'Building' element
        # FIXME assigns all Storeys to all Buildings
        storey_id = 1
        for ifc_building_storey in ifc_file.by_type("IfcBuildingStorey"):
            building_storey = root.createElement("BuildingStorey")
            building_storey.setAttribute(
                "id", fix_xml_stry(ifc_building_storey.GlobalId)
            )
            building.appendChild(building_storey)

            dict_id[fix_xml_stry(ifc_building_storey.GlobalId)] = building_storey

            name = root.createElement("Name")
            name.appendChild(root.createTextNode("Storey_%d" % storey_id))
            storey_id += 1
            building_storey.appendChild(name)

            level = root.createElement("Level")
            level.appendChild(root.createTextNode(str(ifc_building_storey.Elevation)))
            building_storey.appendChild(level)

        # FIXME assigns all Spaces to all Buildings
        # Specify the 'Space' element of the gbXML schema; making use of IFC entity 'IfcSpace'
        # This new element is added as child to the earlier created 'Building' element
        space_id = 1
        for ifc_space in ifc_file.by_type("IfcSpace"):
            space = root.createElement("Space")
            space.setAttribute("id", fix_xml_spc(ifc_space.GlobalId))
            building.appendChild(space)

            dict_id[fix_xml_spc(ifc_space.GlobalId)] = space

            # Refer to the relating 'BuildingStorey' GUID by iterating through IFC entities
            space.setAttribute(
                "buildingStoreyIdRef",
                fix_xml_stry(ifc_space.Decomposes[0].RelatingObject.GlobalId),
            )

            area = root.createElement("Area")
            volume = root.createElement("Volume")

            for ifc_rel in ifc_space.IsDefinedBy:
                if ifc_rel.is_a("IfcRelDefinesByProperties"):
                    if ifc_rel.RelatingPropertyDefinition.is_a("IfcPropertySet"):
                        for (
                            ifc_property
                        ) in ifc_rel.RelatingPropertyDefinition.HasProperties:
                            ifc_value = ifc_property.NominalValue.wrappedValue
                            if ifc_property.Name == "Area":
                                area.appendChild(root.createTextNode(str(ifc_value)))
                                space.appendChild(area)
                            if ifc_property.Name == "Volume":
                                volume.appendChild(root.createTextNode(str(ifc_value)))
                                space.appendChild(volume)

            name = root.createElement("Name")
            name.appendChild(root.createTextNode("Space_%d" % space_id))
            space_id += 1
            space.appendChild(name)

            # Specify the 'SpaceBoundary' element of the gbXML schema; making use of IFC entity 'IfcSpace'
            # This new element is added as child to the earlier created 'Space' element
            for ifc_rel_space_boundary in ifc_space.BoundedBy:

                # Make sure a 'SpaceBoundary' is representing an actual element
                if ifc_rel_space_boundary.RelatedBuildingElement == None:
                    continue

                vertices = get_boundary_vertices(ifc_rel_space_boundary)

                # Create 'SpaceBoundary' elements for the following building elements
                if (
                    ifc_rel_space_boundary.RelatedBuildingElement.is_a("IfcCovering")
                    or ifc_rel_space_boundary.RelatedBuildingElement.is_a("IfcSlab")
                    or ifc_rel_space_boundary.RelatedBuildingElement.is_a("IfcWall")
                    or ifc_rel_space_boundary.RelatedBuildingElement.is_a("IfcRoof")
                ):

                    space_boundary = root.createElement("SpaceBoundary")
                    space_boundary.setAttribute("isSecondLevelBoundary", "true")

                    # Refer to the relating 'SpaceBoundary' GUID by iterating through IFC entities
                    space_boundary.setAttribute(
                        "surfaceIdRef", fix_xml_id(ifc_rel_space_boundary.GlobalId)
                    )

                    space.appendChild(space_boundary)

                    planar_geometry = root.createElement("PlanarGeometry")
                    space_boundary.appendChild(planar_geometry)

                    ifc_relating_space = ifc_rel_space_boundary.RelatingSpace
                    planar_geometry.appendChild(get_poly_loop(root, vertices, ifc_relating_space))

# Specify the 'Surface' element of the gbXML schema; making use of IFC entity 'IfcRelSpaceBoundary'
# This new element is added as child to the earlier created 'Campus' element
opening_id = 1
for ifc_rel_space_boundary in ifc_file.by_type("IfcRelSpaceBoundary"):

    # Make sure a 'SpaceBoundary' is representing an actual element
    if ifc_rel_space_boundary.RelatedBuildingElement == None:
        continue
    if ifc_rel_space_boundary.ConnectionGeometry.SurfaceOnRelatingElement == None:
        continue

    # Specify the 'IfcCurveBoundedPlane' entity which represents the geometry
    vertices = get_boundary_vertices(ifc_rel_space_boundary)

    # Specify each 'Surface' element and set 'SurfaceType' attributes
    if (
        ifc_rel_space_boundary.RelatedBuildingElement.is_a("IfcCovering")
        or ifc_rel_space_boundary.RelatedBuildingElement.is_a("IfcSlab")
        or ifc_rel_space_boundary.RelatedBuildingElement.is_a("IfcWall")
        or ifc_rel_space_boundary.RelatedBuildingElement.is_a("IfcRoof")
    ):

        surface = root.createElement("Surface")
        surface.setAttribute("id", fix_xml_id(ifc_rel_space_boundary.GlobalId))
        dict_id[fix_xml_id(ifc_rel_space_boundary.GlobalId)] = surface

        if ifc_rel_space_boundary.RelatedBuildingElement.is_a("IfcCovering"):
            surface.setAttribute("surfaceType", "Ceiling")

        if ifc_rel_space_boundary.RelatedBuildingElement.is_a("IfcRoof"):
            surface.setAttribute("surfaceType", "Roof")

        if ifc_rel_space_boundary.RelatedBuildingElement.is_a("IfcSlab"):
            if ifc_rel_space_boundary.InternalOrExternalBoundary == "EXTERNAL_EARTH":
                surface.setAttribute("surfaceType", "SlabOnGrade")
            else:
                surface.setAttribute("surfaceType", "InteriorFloor")

        if ifc_rel_space_boundary.RelatedBuildingElement.is_a("IfcWall"):
            if ifc_rel_space_boundary.InternalOrExternalBoundary == "EXTERNAL":
                surface.setAttribute("surfaceType", "ExteriorWall")
            elif ifc_rel_space_boundary.InternalOrExternalBoundary == "EXTERNAL_FIRE":
                surface.setAttribute("surfaceType", "ExteriorWall")
            elif ifc_rel_space_boundary.InternalOrExternalBoundary == "INTERNAL":
                surface.setAttribute("surfaceType", "InteriorWall")
            else:
                surface.setAttribute("surfaceType", "InteriorWall")

        # Refer to the relating 'IfcRelAssociatesMaterial' GUID by iterating through IFC entities
        surface.setAttribute(
            "constructionIdRef",
            fix_xml_cons(
                ifc_rel_space_boundary.RelatedBuildingElement.HasAssociations[
                    0
                ].GlobalId
            ),
        )

        name = root.createElement("Name")
        name.appendChild(
            root.createTextNode(fix_xml_name(ifc_rel_space_boundary.GlobalId))
        )

        surface.appendChild(name)

        adjacent_space_id = root.createElement("AdjacentSpaceId")

        # Refer to the relating 'Space' GUID by iterating through IFC entities
        adjacent_space_id.setAttribute(
            "spaceIdRef", fix_xml_spc(ifc_rel_space_boundary.RelatingSpace.GlobalId)
        )
        surface.appendChild(adjacent_space_id)

        planar_geometry = root.createElement("PlanarGeometry")
        surface.appendChild(planar_geometry)

        ifc_relating_space = ifc_rel_space_boundary.RelatingSpace
        planar_geometry.appendChild(get_poly_loop(root, vertices, ifc_relating_space))

        cad_object_id = root.createElement("CADObjectId")
        cad_object_id.appendChild(
            root.createTextNode(fix_xml_name(ifc_rel_space_boundary.GlobalId))
        )
        surface.appendChild(cad_object_id)

        campus.appendChild(surface)

    if ifc_rel_space_boundary.RelatedBuildingElement.is_a(
        "IfcWindow"
    ) or ifc_rel_space_boundary.RelatedBuildingElement.is_a("IfcDoor"):
        opening = root.createElement("Opening")

        # Refer to the relating 'IfcWindow' GUID by iterating through IFC entities
        opening.setAttribute(
            "windowTypeIdRef",
            fix_xml_id(ifc_rel_space_boundary.RelatedBuildingElement.GlobalId),
        )
        if ifc_rel_space_boundary.RelatedBuildingElement.is_a("IfcWindow"):
            opening.setAttribute("openingType", "OperableWindow")
        elif ifc_rel_space_boundary.RelatedBuildingElement.is_a("IfcDoor"):
            opening.setAttribute("openingType", "NonSlidingDoor")

        opening.setAttribute("id", "Opening%d" % opening_id)
        opening_id += 1

        # If the building element is an 'IfcWindow' the gbXML element 'Opening' is added
        planar_geometry = root.createElement("PlanarGeometry")
        opening.appendChild(planar_geometry)

        ifc_relating_space = ifc_rel_space_boundary.RelatingSpace
        new_z = ifc_relating_space.ObjectPlacement.PlacementRelTo.RelativePlacement.Location.Coordinates[
            2
        ]
        planar_geometry.appendChild(get_poly_loop(root, vertices, new_z))

        name = root.createElement("Name")
        name.appendChild(
            root.createTextNode(
                fix_xml_name(ifc_rel_space_boundary.RelatedBuildingElement.Name)
            )
        )
        opening.appendChild(name)

        cad_object_id = root.createElement("CADObjectId")
        cad_object_id.appendChild(
            root.createTextNode(
                fix_xml_name(ifc_rel_space_boundary.RelatedBuildingElement.Name)
            )
        )
        opening.appendChild(cad_object_id)

        surface.appendChild(opening)

    else:
        continue

# Specify the 'WindowType' element of the gbXML schema; making use of IFC entity 'IfcWindow'
# This new element is added as child to the earlier created 'gbXML' element
for ifc_building_element in ifc_file.by_type("IfcWindow") + ifc_file.by_type("IfcDoor"):

    # There doesn't appear to be a 'DoorType'
    if ifc_building_element.is_a("IfcWindow") or ifc_building_element.is_a("IfcDoor"):
        window_type = root.createElement("WindowType")

    window_type.setAttribute("id", fix_xml_id(ifc_building_element.GlobalId))
    gbxml.appendChild(window_type)

    dict_id[fix_xml_id(ifc_building_element.GlobalId)] = window_type

    name = root.createElement("Name")
    name.appendChild(root.createTextNode(fix_xml_name(ifc_building_element.Name)))
    window_type.appendChild(name)

    description = root.createElement("Description")
    description.appendChild(
        root.createTextNode(fix_xml_name(ifc_building_element.Name))
    )
    window_type.appendChild(description)

    # Specify analytical properties of the 'IfcWindow' by iterating through IFC entities

    u_value = root.createElement("U-value")
    for ifc_rel in ifc_building_element.IsDefinedBy:
        if ifc_rel.is_a("IfcRelDefinesByProperties"):
            if ifc_rel.RelatingPropertyDefinition.is_a("IfcPropertySet"):
                for ifc_property in ifc_rel.RelatingPropertyDefinition.HasProperties:
                    if ifc_property.Name == "ThermalTransmittance":
                        ifc_value = ifc_property.NominalValue.wrappedValue
                        u_value.setAttribute("unit", "WPerSquareMeterK")
                        u_value.appendChild(root.createTextNode(str(ifc_value)))
                        window_type.appendChild(u_value)

    solar_heat_gain_coeff = root.createElement("SolarHeatGainCoeff")
    transmittance = root.createElement("Transmittance")
    for ifc_rel in ifc_building_element.IsDefinedBy:
        if ifc_rel.is_a("IfcRelDefinesByType"):
            if ifc_rel.RelatingType.is_a("IfcTypeProduct"):
                for ifc_property_set in ifc_rel.RelatingType.HasPropertySets:
                    if ifc_property_set.Name == "Analytical Properties(Type)":
                        for ifc_property in ifc_property_set.HasProperties:
                            ifc_value = ifc_property.NominalValue.wrappedValue
                            if ifc_property.Name == "Solar Heat Gain Coefficient":
                                solar_heat_gain_coeff.setAttribute("unit", "Fraction")
                                solar_heat_gain_coeff.appendChild(
                                    root.createTextNode(str(ifc_value))
                                )
                                window_type.appendChild(solar_heat_gain_coeff)

                            if ifc_property.Name == "Visual Light Transmittance":
                                transmittance.setAttribute("unit", "Fraction")
                                transmittance.setAttribute("type", "Visible")
                                transmittance.appendChild(
                                    root.createTextNode(str(ifc_value))
                                )
                                window_type.appendChild(transmittance)

# Specify the 'Construction' element of the gbXML schema; making use of IFC entity 'IfcRelSpaceBoundary'
# This new element is added as child to the earlier created 'gbXML' element
ifc_global_ids = []

for ifc_rel_space_boundary in ifc_file.by_type("IfcRelSpaceBoundary"):
    # Make sure a 'SpaceBoundary' is representing an actual element
    if ifc_rel_space_boundary.RelatedBuildingElement is None:
        continue

    if (
        ifc_rel_space_boundary.RelatedBuildingElement.is_a("IfcCovering")
        or ifc_rel_space_boundary.RelatedBuildingElement.is_a("IfcSlab")
        or ifc_rel_space_boundary.RelatedBuildingElement.is_a("IfcWall")
        or ifc_rel_space_boundary.RelatedBuildingElement.is_a("IfcRoof")
    ):

        # Refer to the relating 'IfcRelAssociatesMaterial' GUID by iterating through IFC entities
        ifc_global_id = ifc_rel_space_boundary.RelatedBuildingElement.HasAssociations[
            0
        ].GlobalId

        # Make use of a list to make sure no same 'Construction' elements are added twice
        if ifc_global_id not in ifc_global_ids:
            ifc_global_ids.append(ifc_global_id)

            construction = root.createElement("Construction")
            construction.setAttribute("id", fix_xml_cons(ifc_global_id))
            dict_id[fix_xml_cons(ifc_global_id)] = construction

            # Specify analytical properties of the 'Construction' element by iterating through IFC entities

            u_value = root.createElement("U-value")
            for ifc_rel in ifc_rel_space_boundary.RelatedBuildingElement.IsDefinedBy:
                if ifc_rel.is_a("IfcRelDefinesByProperties"):
                    if ifc_rel.RelatingPropertyDefinition.is_a("IfcPropertySet"):
                        for (
                            ifc_property
                        ) in ifc_rel.RelatingPropertyDefinition.HasProperties:
                            ifc_value = ifc_property.NominalValue.wrappedValue

                            if ifc_rel_space_boundary.RelatedBuildingElement.is_a(
                                "IfcWall"
                            ):
                                if ifc_property.Name == "ThermalTransmittance":
                                    u_value.setAttribute("unit", "WPerSquareMeterK")
                                    u_value.appendChild(
                                        root.createTextNode(str(ifc_value))
                                    )
                                    construction.appendChild(u_value)

                            if ifc_property.Name == "Heat Transfer Coefficient (U)":
                                u_value.setAttribute("unit", "WPerSquareMeterK")
                                u_value.appendChild(root.createTextNode(str(ifc_value)))
                                construction.appendChild(u_value)

            absorptance = root.createElement("Absorptance")
            for ifc_rel in ifc_rel_space_boundary.RelatedBuildingElement.IsDefinedBy:
                if ifc_rel.is_a("IfcRelDefinesByProperties"):
                    if ifc_rel.RelatingPropertyDefinition.is_a("IfcPropertySet"):
                        for (
                            ifc_property
                        ) in ifc_rel.RelatingPropertyDefinition.HasProperties:
                            if ifc_property.Name == "Absorptance":
                                ifc_value = ifc_property.NominalValue.wrappedValue
                                absorptance.setAttribute("unit", "Fraction")
                                absorptance.setAttribute("type", "ExtIR")
                                absorptance.appendChild(
                                    root.createTextNode(str(ifc_value))
                                )
                                construction.appendChild(absorptance)

            for (
                association
            ) in ifc_rel_space_boundary.RelatedBuildingElement.HasAssociations:
                if association.is_a("IfcRelAssociatesMaterial"):
                    # Refer to the relating 'IfcRelAssociatesMaterial' GUID by iterating through IFC entities

                    layer_id = root.createElement("LayerId")
                    layer_id.setAttribute(
                        "layerIdRef", fix_xml_layer(association.GlobalId)
                    )
                    construction.appendChild(layer_id)

                    # Refer to the relating 'IfcMaterialLayerSet' name by iterating through IFC entities
                    name = root.createElement("Name")
                    if hasattr(association.RelatingMaterial, "ForLayerSet"):
                        name.appendChild(
                            root.createTextNode(
                                association.RelatingMaterial.ForLayerSet.LayerSetName
                            )
                        )
                    construction.appendChild(name)

            gbxml.appendChild(construction)

    else:
        continue

# Specify the 'Layer' element of the gbXML schema; making use of IFC entity 'IfcBuildingElement'
# This new element is added as child to the earlier created 'gbXML' element
ifc_building_elements = ifc_file.by_type("IfcBuildingElement")
for ifc_building_element in ifc_building_elements:
    if (
        ifc_building_element.is_a("IfcWall")
        or ifc_building_element.is_a("IfcCovering")
        or ifc_building_element.is_a("IfcSlab")
        or ifc_building_element.is_a("IfcRoof")
    ):

        # Try and catch an Element that is just an Aggregate
        if ifc_building_element.IsDecomposedBy:
            continue
        for association in ifc_building_element.HasAssociations:
            if association.is_a("IfcRelAssociatesMaterial"):
                # Refer to the relating 'IfcRelAssociatesMaterial' GUID by iterating through IFC entities

                layer = root.createElement("Layer")
                layer.setAttribute("id", fix_xml_layer(association.GlobalId))

                dict_id[fix_xml_layer(association.GlobalId)] = layer

                # Specify the 'IfcMaterialLayer' entity and iterate to each 'IfcMaterial' entity
                if not association.RelatingMaterial.is_a("IfcMaterialLayerSetUsage"):
                    continue
                for (
                    ifc_material_layer
                ) in association.RelatingMaterial.ForLayerSet.MaterialLayers:
                    material_id = root.createElement("MaterialId")
                    material_id.setAttribute(
                        "materialIdRef", "mat_%d" % ifc_material_layer.Material.id()
                    )
                    layer.appendChild(material_id)

                    dict_id["mat_%d" % ifc_material_layer.Material.id()] = layer

                gbxml.appendChild(layer)

    else:
        continue

# Specify the 'Material' element of the gbXML schema; making use of IFC entity 'IfcBuildingElement'
# This new element is added as child to the earlier created 'gbXML' element
ifc_material_layer_ids = []

for ifc_building_element in ifc_building_elements:
    if (
        ifc_building_element.is_a("IfcWall")
        or ifc_building_element.is_a("IfcSlab")
        or ifc_building_element.is_a("IfcCovering")
        or ifc_building_element.is_a("IfcRoof")
    ):

        # Try and catch an Element that is just an Aggregate
        if ifc_building_element.IsDecomposedBy:
            continue

        for association in ifc_building_element.HasAssociations:
            if not association.is_a("IfcRelAssociatesMaterial"):
                continue
            if not association.RelatingMaterial.is_a("IfcMaterialLayerSetUsage"):
                continue

            for (
                ifc_material_layer
            ) in association.RelatingMaterial.ForLayerSet.MaterialLayers:
                ifc_material_layer_id = ifc_material_layer.Material.id()

                # Make use of a list to make sure no same 'Materials' elements are added twice
                if ifc_material_layer_id not in ifc_material_layer_ids:
                    ifc_material_layer_ids.append(ifc_material_layer_id)

                    material = root.createElement("Material")
                    material.setAttribute("id", "mat_%d" % ifc_material_layer_id)
                    dict_id["mat_%d" % ifc_material_layer_id] = material

                    name = root.createElement("Name")
                    name.appendChild(
                        root.createTextNode(ifc_material_layer.Material.Name)
                    )
                    material.appendChild(name)

                    thickness = root.createElement("Thickness")
                    thickness.setAttribute("unit", "Meters")
                    layer_thickness = ifc_material_layer.LayerThickness
                    thickness.appendChild(root.createTextNode((str(layer_thickness))))
                    material.appendChild(thickness)

                    r_value = root.createElement("R-value")
                    r_value.setAttribute("unit", "SquareMeterKPerW")

                    # Analytical properties of the Material entity can be found directly
                    if hasattr(ifc_material_layer.Material, "HasProperties"):
                        for (
                            ifc_material_property
                        ) in ifc_material_layer.Material.HasProperties:
                            if ifc_material_property.Name == "Pset_MaterialEnergy":
                                for (
                                    pset_material_energy
                                ) in ifc_material_property.Properties:
                                    if (
                                        pset_material_energy.Name
                                        == "ThermalConductivityTemperatureDerivative"
                                    ):
                                        ifc_value = (
                                            pset_material_energy.NominalValue.wrappedValue
                                        )
                                        r_value.setAttribute("unit", "SquareMeterKPerW")
                                        r_value.appendChild(
                                            root.createTextNode(str(ifc_value))
                                        )
                                        material.appendChild(r_value)

                                        gbxml.appendChild(material)

                    # Specify analytical properties of the 'Material' element by iterating through IFC entities
                    for ifc_rel in ifc_building_element.IsDefinedBy:
                        if ifc_rel.is_a("IfcRelDefinesByType"):
                            if ifc_rel.RelatingType.is_a("IfcWallType"):
                                for (
                                    ifc_property_set
                                ) in ifc_rel.RelatingType.HasPropertySets:
                                    if (
                                        ifc_property_set.Name
                                        == "Analytical Properties(Type)"
                                    ):
                                        for (
                                            ifc_property
                                        ) in ifc_property_set.HasProperties:
                                            if (
                                                ifc_property.Name
                                                == "Heat Transfer Coefficient (U)"
                                            ):
                                                ifc_value = (
                                                    ifc_property.NominalValue.wrappedValue
                                                )
                                                r_value.appendChild(
                                                    root.createTextNode(
                                                        str(layer_thickness / ifc_value)
                                                    )
                                                )
                                                material.appendChild(r_value)

                                                gbxml.appendChild(material)

                        elif ifc_rel.is_a("IfcRelDefinesByProperties"):
                            if ifc_rel.RelatingPropertyDefinition.is_a(
                                "IfcPropertySet"
                            ):
                                for (
                                    ifc_property
                                ) in ifc_rel.RelatingPropertyDefinition.HasProperties:
                                    if (
                                        ifc_property.Name
                                        == "Heat Transfer Coefficient (U)"
                                    ):
                                        ifc_value = (
                                            ifc_property.NominalValue.wrappedValue
                                        )
                                        r_value.setAttribute("unit", "SquareMeterKPerW")
                                        r_value.appendChild(
                                            root.createTextNode(
                                                str(layer_thickness / ifc_value)
                                            )
                                        )
                                        material.appendChild(r_value)

                                        gbxml.appendChild(material)

                        if ifc_building_element.is_a("IfcCovering"):
                            if ifc_rel.is_a("IfcRelDefinesByProperties") and hasattr(
                                ifc_rel, "RelatingType"
                            ):
                                if ifc_rel.RelatingType.is_a("IfcPropertySet"):
                                    for (
                                        ifc_property_set
                                    ) in ifc_rel.RelatingType.HasPropertySets:
                                        if (
                                            ifc_property_set.Name
                                            == "Analytical Properties(Type)"
                                        ):
                                            for (
                                                ifc_property
                                            ) in ifc_property_set.HasProperties:
                                                if (
                                                    ifc_property.Name
                                                    == "Heat Transfer Coefficient (U)"
                                                ):
                                                    ifc_value = (
                                                        ifc_property.NominalValue.wrappedValue
                                                    )
                                                    r_value.setAttribute(
                                                        "unit", "SquareMeterKPerW"
                                                    )
                                                    r_value.appendChild(
                                                        root.createTextNode(
                                                            str(
                                                                layer_thickness
                                                                / ifc_value
                                                            )
                                                        )
                                                    )
                                                    material.appendChild(r_value)

                                                    gbxml.appendChild(material)


def create_DocumentHistory(ifc_file, root):
    """Specify the 'DocumentHistory' element of the gbXML schema; making use of IFC entity 'IfcApplication' and 'IfcPerson'"""
    document_history = root.createElement("DocumentHistory")
    for ifc_application in ifc_file.by_type("IfcApplication"):

        program_info = root.createElement("ProgramInfo")
        program_info.setAttribute("id", ifc_application.ApplicationIdentifier)

        company_name = root.createElement("CompanyName")
        company_name.appendChild(
            root.createTextNode(ifc_application.ApplicationDeveloper.Name)
        )
        program_info.appendChild(company_name)

        product_name = root.createElement("ProductName")
        product_name.appendChild(
            root.createTextNode(ifc_application.ApplicationFullName)
        )
        program_info.appendChild(product_name)

        version = root.createElement("Version")
        version.appendChild(root.createTextNode(ifc_application.Version))
        program_info.appendChild(version)

        document_history.appendChild(program_info)

    for ifc_person in ifc_file.by_type("IfcPerson"):

        created_by = root.createElement("CreatedBy")
        created_by.setAttribute("personId", ifc_person.GivenName)
        for ifc_application in ifc_file.by_type("IfcApplication"):
            created_by.setAttribute("programId", ifc_application.ApplicationIdentifier)
        today = datetime.date.today()
        created_by.setAttribute(
            "date", today.strftime("%Y-%m-%dT") + time.strftime("%H:%M:%S")
        )

        document_history.appendChild(created_by)

        person_info = root.createElement("PersonInfo")
        person_info.setAttribute("id", ifc_person.GivenName)

        document_history.appendChild(person_info)

    return document_history


gbxml.appendChild(create_DocumentHistory(ifc_file, root))

# Create a new XML file and write all created elements to it
save_path_file = sys.argv[2]

root.writexml(open(save_path_file, "w"), indent="  ", addindent="  ", newl="\n")
