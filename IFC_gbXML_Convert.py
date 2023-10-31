#!/usr/bin/python3

import ifcopenshell.util.placement
import numpy as np
import datetime
import time
import sys
from xml.dom import minidom

get_pset = ifcopenshell.util.element.get_pset

# Copyright (C) 2016-2023
# Maarten Visschers <maartenvisschers@hotmail.com>
# Bruno Postle <bruno@postle.net>
# Tokarzewski <bartlomiej.tokarzewski@gmail.com>
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# this program. If not, see <https://www.gnu.org/licenses/>.

# from https://github.com/wassimj/TopologicSverchok/blob/main/nodes/Topologic/ifc_topologic.py
def get_boundary_vertices(ifc_rel_space_boundary):
    ifc_curve = ifc_rel_space_boundary.ConnectionGeometry.SurfaceOnRelatingElement
    if ifc_curve.is_a("IfcFaceSurface"):
        # FIXME assumes Bound is IfcPolyLoop, could also be IfcEdgeLoop or IfcVertexLoop
        ifc_points = ifc_curve.Bounds[0].Bound.Polygon
        plane_coords = [v.Coordinates for v in ifc_points]
        ifc_plane = ifc_curve.FaceSurface
    elif ifc_curve.is_a("IfcCurveBoundedPlane"):
        # assumes OuterBoundary is IfcIndexedPolyCurve or IfcPolyline
        if ifc_curve.OuterBoundary.is_a("IfcPolyline"):
            plane_coords = [
                point.Coordinates for point in ifc_curve.OuterBoundary.Points
            ]
        elif ifc_curve.OuterBoundary.is_a("IfcIndexedPolyCurve"):
            plane_coords = ifc_curve.OuterBoundary.Points[0]
        else:
            return []
        ifc_plane = ifc_curve.BasisSurface
    space_matrix = ifcopenshell.util.placement.get_local_placement(
        ifc_rel_space_boundary.RelatingSpace.ObjectPlacement
    )
    plane_matrix = ifcopenshell.util.placement.get_axis2placement(ifc_plane.Position)
    plane_vertices = [
        np.array(v + (0, 1)) if len(v) == 2 else np.array(v + (1,))
        for v in plane_coords
    ]
    vertices = [(space_matrix @ plane_matrix @ v)[0:3] for v in plane_vertices]
    if (
        np.dot(
            np.cross(vertices[1] - vertices[0], vertices[-1] - vertices[0]),
            (plane_matrix @ np.array([0, 0, 1, 0]))[0:3],
        )
        < -1e-6
    ):
        vertices.reverse()
    return vertices


def get_poly_loop(root, vertices):
    poly_loop = root.createElement("PolyLoop")
    # gbXML PolyLoops don't have coincident start/end vertices
    if np.allclose(vertices[0], vertices[-1]):
        del vertices[-1]
    for v in vertices:
        x, y, z = v
        cartesian_point = root.createElement("CartesianPoint")
        for c in x, y, z:
            if "e-" in str(c):
                c = 0.0
            coordinate = root.createElement("Coordinate")
            coordinate.appendChild(root.createTextNode(str(c)))
            cartesian_point.appendChild(coordinate)
        poly_loop.appendChild(cartesian_point)
    return poly_loop


# valid GUID chars: ^[0123][0-9a-zA-Z_$]{21}$
# valid XML id chars: ^[a-zA-Z_][a-zA-Z0-9_.-]$

# Align the gbXML input according to the predefined official gbXML schema
def remove_unnecessary_characters(element):
    char_to_replace = {"$": "", ":": "", " ": "", "(": "", ")": ""}
    for key, value in char_to_replace.items():
        element = element.replace(key, value)
    return element


def fix_xml_cmps(element):
    return "campus_" + element.replace("$", "-")


def fix_xml_bldng(element):
    return "building_" + element.replace("$", "-")


def fix_xml_stry(element):
    return "buildingstorey_" + element.replace("$", "-")


def fix_xml_spc(element):
    return "space_" + element.replace("$", "-")


def fix_xml_id(element):
    return "id_" + element.replace("$", "-")


# Only fix_xml_name() is used for non-guid data
def fix_xml_name(element):
    return "object_" + remove_unnecessary_characters(element)


def fix_xml_cons(element):
    return "construction_" + element.replace("$", "-")


def fix_xml_layer(element):
    return "layer_" + element.replace("$", "-")


def create_gbxml(ifc_file):
    # Create the XML root by making use of MiniDom
    root = minidom.Document()

    # Create the 'gbXML' element and append it to the Root of the document
    gbxml = root.createElement("gbXML")
    root.appendChild(gbxml)

    # FIXME support non-SI units
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

    for ifc_site in ifc_file.by_type("IfcSite"):
        if ifc_site.RefLongitude == None:
            ifc_site.RefLatitude = [53, 23, 0]
            ifc_site.RefLongitude = [1, 28, 0]
            ifc_site.RefElevation = 75.0
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

        ifc_postal_address = ifc_site.SiteAddress
        zipcode = root.createElement("ZipcodeOrPostalCode")
        location.appendChild(zipcode)
        if ifc_postal_address:
            zipcode.appendChild(root.createTextNode(ifc_postal_address.PostalCode))

            name = root.createElement("Name")
            name.appendChild(
                root.createTextNode(
                    str(ifc_postal_address.Region)
                    + ", "
                    + str(ifc_postal_address.Country)
                )
            )
            location.appendChild(name)
        else:
            zipcode.appendChild(root.createTextNode("S1 2GH"))

        # Specify the 'Building' element of the gbXML schema; making use of IFC entity 'IfcBuilding'
        # This new element is added as child to the earlier created 'Campus' element
        if not ifc_site.IsDecomposedBy:
            continue

        for ifc_building in ifc_site.IsDecomposedBy[0].RelatedObjects:
            if not ifc_building.is_a("IfcBuilding"):
                continue
            building = root.createElement("Building")
            building.setAttribute("id", fix_xml_bldng(ifc_building.GlobalId))
            building.setAttribute("buildingType", "Unknown")
            campus.appendChild(building)

            dict_id[fix_xml_bldng(ifc_building.GlobalId)] = building

            ifc_postal_address = ifc_building.BuildingAddress
            if ifc_postal_address:
                street_address = root.createElement("StreetAddress")
                street_address.appendChild(
                    root.createTextNode(
                        str(ifc_postal_address.Region)
                        + ", "
                        + str(ifc_postal_address.Country)
                    )
                )
                building.appendChild(street_address)

            # Specify the 'BuildingStorey' element of the gbXML schema; making use of IFC entity 'IfcBuildingStorey'
            # This new element is added as child to the earlier created 'Building' element
            if not ifc_building.IsDecomposedBy:
                continue

            for ifc_building_storey in ifc_building.IsDecomposedBy[0].RelatedObjects:
                if not ifc_building_storey.is_a("IfcBuildingStorey"):
                    continue
                building_storey = root.createElement("BuildingStorey")
                building_storey.setAttribute(
                    "id", fix_xml_stry(ifc_building_storey.GlobalId)
                )
                building.appendChild(building_storey)

                dict_id[fix_xml_stry(ifc_building_storey.GlobalId)] = building_storey

                name = root.createElement("Name")
                name.appendChild(root.createTextNode(ifc_building_storey.Name))
                building_storey.appendChild(name)

                level = root.createElement("Level")
                level.appendChild(
                    root.createTextNode(
                        str(
                            ifcopenshell.util.placement.get_storey_elevation(
                                ifc_building_storey
                            )
                        )
                    )
                )
                building_storey.appendChild(level)

                # Specify the 'Space' element of the gbXML schema; making use of IFC entity 'IfcSpace'
                # This new element is added as child to the earlier created 'Building' element
                if not ifc_building_storey.IsDecomposedBy:
                    continue

                for ifc_space in ifc_building_storey.IsDecomposedBy[0].RelatedObjects:
                    if not ifc_space.is_a("IfcSpace"):
                        continue
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
                    # area.setAttribute("unit", "MetreSquare")
                    area.appendChild(root.createTextNode("1.0"))
                    space.appendChild(area)
                    volume = root.createElement("Volume")
                    # volume.setAttribute("unit", "MetreCube)
                    volume.appendChild(root.createTextNode("1.0"))
                    space.appendChild(volume)

                    pset_area = (
                        get_pset(
                            # IFC4
                            ifc_space,
                            "Qto_SpaceBaseQuantities",
                            prop="NetFloorArea",
                        )
                        or get_pset(
                            # IFC4
                            ifc_space,
                            "Pset_SpaceCommon",
                            prop="NetPlannedArea",
                        )
                        or get_pset(
                            # IFC2X3
                            ifc_space,
                            "BaseQuantities",
                            prop="NetFloorArea",
                        )
                        or get_pset(
                            # IFC2X3
                            ifc_space,
                            "Dimensions",
                            prop="Area",
                        )
                    )
                    if pset_area:
                        area.firstChild.data = str(pset_area)

                    pset_volume = (
                        get_pset(
                            # IFC4
                            ifc_space,
                            "Qto_SpaceBaseQuantities",
                            prop="NetVolume",
                        )
                        or get_pset(
                            # IFC2X3
                            ifc_space,
                            "BaseQuantities",
                            prop="GrossVolume",
                        )
                        or get_pset(
                            # IFC2X3
                            ifc_space,
                            "Dimensions",
                            prop="Volume",
                        )
                    )
                    if pset_volume:
                        volume.firstChild.data = str(pset_volume)

                    name = root.createElement("Name")
                    name.appendChild(root.createTextNode(ifc_space.Name))
                    space.appendChild(name)

                    # Specify the 'SpaceBoundary' element of the gbXML schema; making use of IFC entity 'IfcSpace'
                    # This new element is added as child to the earlier created 'Space' element
                    for ifc_rel_space_boundary in ifc_space.BoundedBy:

                        # Make sure a 'SpaceBoundary' is representing an actual element
                        if ifc_rel_space_boundary.RelatedBuildingElement == None:
                            continue

                        vertices = get_boundary_vertices(ifc_rel_space_boundary)

                        # Create 'SpaceBoundary' elements for the following building elements
                        if ifc_rel_space_boundary.RelatedBuildingElement.is_a() in [
                            "IfcWall",
                            "IfcSlab",
                            "IfcRoof",
                            "IfcCovering",
                        ]:

                            space_boundary = root.createElement("SpaceBoundary")
                            space_boundary.setAttribute("isSecondLevelBoundary", "true")

                            # Refer to the relating 'SpaceBoundary' GUID by iterating through IFC entities
                            space_boundary.setAttribute(
                                "surfaceIdRef",
                                fix_xml_id(ifc_rel_space_boundary.GlobalId),
                            )

                            space.appendChild(space_boundary)

                            planar_geometry = root.createElement("PlanarGeometry")
                            space_boundary.appendChild(planar_geometry)

                            planar_geometry.appendChild(get_poly_loop(root, vertices))

    # Specify the 'Surface' element of the gbXML schema; making use of IFC entity 'IfcRelSpaceBoundary'
    # This new element is added as child to the earlier created 'Campus' element
    for ifc_rel_space_boundary in ifc_file.by_type("IfcRelSpaceBoundary"):

        # Make sure a 'SpaceBoundary' is representing an actual element
        if ifc_rel_space_boundary.RelatedBuildingElement == None:
            continue
        if ifc_rel_space_boundary.ConnectionGeometry.SurfaceOnRelatingElement == None:
            continue
        if ifc_rel_space_boundary.RelatingSpace == None:
            continue

        # Specify the 'IfcCurveBoundedPlane' entity which represents the geometry
        vertices = get_boundary_vertices(ifc_rel_space_boundary)

        # Specify each 'Surface' element and set 'SurfaceType' attributes
        if ifc_rel_space_boundary.RelatedBuildingElement.is_a() in [
            "IfcWall",
            "IfcSlab",
            "IfcRoof",
            "IfcCovering",
        ]:

            surface = root.createElement("Surface")
            surface.setAttribute("id", fix_xml_id(ifc_rel_space_boundary.GlobalId))
            dict_id[fix_xml_id(ifc_rel_space_boundary.GlobalId)] = surface

            if ifc_rel_space_boundary.RelatedBuildingElement.is_a("IfcCovering"):
                surface.setAttribute("surfaceType", "Ceiling")

            if ifc_rel_space_boundary.RelatedBuildingElement.is_a("IfcRoof"):
                surface.setAttribute("surfaceType", "Roof")
                surface.setAttribute("exposedToSun", "true")

            if ifc_rel_space_boundary.RelatedBuildingElement.is_a("IfcSlab"):
                if (
                    ifc_rel_space_boundary.InternalOrExternalBoundary
                    == "EXTERNAL_EARTH"
                ):
                    surface.setAttribute("surfaceType", "SlabOnGrade")
                else:
                    surface.setAttribute("surfaceType", "InteriorFloor")

            if ifc_rel_space_boundary.RelatedBuildingElement.is_a("IfcWall"):
                if ifc_rel_space_boundary.InternalOrExternalBoundary == "EXTERNAL":
                    surface.setAttribute("surfaceType", "ExteriorWall")
                    surface.setAttribute("exposedToSun", "true")
                elif (
                    ifc_rel_space_boundary.InternalOrExternalBoundary == "EXTERNAL_FIRE"
                ):
                    surface.setAttribute("surfaceType", "InteriorWall")
                    surface.setAttribute("exposedToSun", "false")
                else:
                    surface.setAttribute("surfaceType", "InteriorWall")
                    surface.setAttribute("exposedToSun", "false")

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

            planar_geometry.appendChild(get_poly_loop(root, vertices))

            cad_object_id = root.createElement("CADObjectId")
            cad_object_id.appendChild(
                root.createTextNode(fix_xml_name(ifc_rel_space_boundary.GlobalId))
            )
            surface.appendChild(cad_object_id)

            campus.appendChild(surface)

    opening_id = 1
    for ifc_rel_space_boundary in ifc_file.by_type("IfcRelSpaceBoundary"):

        ifc_building_element = ifc_rel_space_boundary.RelatedBuildingElement
        # Make sure a 'SpaceBoundary' is representing an actual element
        if ifc_building_element == None:
            continue
        if ifc_rel_space_boundary.ConnectionGeometry.SurfaceOnRelatingElement == None:
            continue
        if ifc_rel_space_boundary.RelatingSpace == None:
            continue

        # Specify the 'IfcCurveBoundedPlane' entity which represents the geometry
        vertices = get_boundary_vertices(ifc_rel_space_boundary)

        if ifc_building_element.is_a() in [
            "IfcWindow",
            "IfcDoor",
        ]:

            if ifc_building_element.IsTypedBy:
                ifc_building_element = ifc_building_element.IsTypedBy[0].RelatingType

            opening = root.createElement("Opening")
            dict_id[fix_xml_id(ifc_rel_space_boundary.GlobalId)] = opening

            # Refer to the relating 'IfcWindow' GUID by iterating through IFC entities
            opening.setAttribute(
                "windowTypeIdRef",
                fix_xml_id(ifc_building_element.GlobalId),
            )
            if ifc_building_element.is_a() in ["IfcWindow", "IfcWindowType"]:
                opening.setAttribute("openingType", "OperableWindow")
            elif ifc_building_element.is_a() in ["IfcDoor", "IfcDoorType"]:
                opening.setAttribute("openingType", "NonSlidingDoor")

            opening.setAttribute("id", "Opening%d" % opening_id)
            opening_id += 1

            # If the building element is an 'IfcWindow' the gbXML element 'Opening' is added
            planar_geometry = root.createElement("PlanarGeometry")
            opening.appendChild(planar_geometry)

            planar_geometry.appendChild(get_poly_loop(root, vertices))

            name = root.createElement("Name")
            name.appendChild(
                root.createTextNode(fix_xml_name(ifc_building_element.Name))
            )
            opening.appendChild(name)

            cad_object_id = root.createElement("CADObjectId")
            cad_object_id.appendChild(
                root.createTextNode(fix_xml_name(ifc_building_element.Name))
            )
            opening.appendChild(cad_object_id)

            if fix_xml_id(ifc_rel_space_boundary.ParentBoundary.GlobalId) in dict_id:
                surface = dict_id[
                    fix_xml_id(ifc_rel_space_boundary.ParentBoundary.GlobalId)
                ]
                surface.appendChild(opening)

    ifc_building_element_guids = []
    # Specify the 'WindowType' element of the gbXML schema; making use of IFC entity 'IfcWindow'
    # This new element is added as child to the earlier created 'gbXML' element
    for ifc_building_element in ifc_file.by_type("IfcWindow") + ifc_file.by_type(
        "IfcDoor"
    ):

        if ifc_building_element.IsTypedBy:
            ifc_building_element = ifc_building_element.IsTypedBy[0].RelatingType

        ifc_building_element_guid = ifc_building_element.GlobalId
        if ifc_building_element_guid not in ifc_building_element_guids:
            ifc_building_element_guids.append(ifc_building_element_guid)
            # There doesn't appear to be a 'DoorType'
            window_type = root.createElement("WindowType")
            window_type.setAttribute("id", fix_xml_id(ifc_building_element.GlobalId))
            dict_id[fix_xml_id(ifc_building_element.GlobalId)] = window_type

            name = root.createElement("Name")
            name.appendChild(
                root.createTextNode(fix_xml_name(ifc_building_element.Name))
            )
            window_type.appendChild(name)

            description = root.createElement("Description")
            description.appendChild(
                root.createTextNode(fix_xml_name(ifc_building_element.Name))
            )
            window_type.appendChild(description)

            # Specify analytical properties of the 'IfcWindow' by iterating through IFC entities

            u_value = root.createElement("U-value")
            # FIXME SI units
            u_value.setAttribute("unit", "WPerSquareMeterK")
            u_value.appendChild(root.createTextNode("10.0"))

            pset_u_value = get_pset(
                # IFC4
                ifc_building_element,
                "Pset_DoorCommon",
                prop="ThermalTransmittance",
            ) or get_pset(
                # IFC4
                ifc_building_element,
                "Pset_WindowCommon",
                prop="ThermalTransmittance",
            )
            if pset_u_value:
                u_value.firstChild.data = str(pset_u_value)
                window_type.appendChild(u_value)

            solar_heat_gain_coeff = root.createElement("SolarHeatGainCoeff")
            solar_heat_gain_coeff.setAttribute("unit", "Fraction")
            solar_heat_gain_coeff.appendChild(root.createTextNode("1.0"))

            pset_solar_heat = get_pset(
                # IFC2X3
                ifc_building_element,
                "Analytical Properties(Type)",
                prop="Solar Heat Gain Coefficient",
            )
            if pset_solar_heat:
                solar_heat_gain_coeff.firstChild.data = str(pset_solar_heat)
                window_type.appendChild(solar_heat_gain_coeff)

            transmittance = root.createElement("Transmittance")
            transmittance.setAttribute("unit", "Fraction")
            transmittance.setAttribute("type", "Visible")
            transmittance.appendChild(root.createTextNode("1.0"))

            pset_transmittance = (
                get_pset(
                    # IFC4
                    ifc_building_element,
                    "Pset_DoorCommon",
                    prop="GlazingAreaFraction",
                )
                or get_pset(
                    # IFC4
                    ifc_building_element,
                    "Pset_WindowCommon",
                    prop="GlazingAreaFraction",
                )
                or get_pset(
                    # IFC2X3
                    ifc_building_element,
                    "Analytical Properties(Type)",
                    prop="Visual Light Transmittance",
                )
            )
            if pset_transmittance:
                transmittance.firstChild.data = str(pset_transmittance)
                window_type.appendChild(transmittance)

            gbxml.appendChild(window_type)

    # Specify the 'Construction' element of the gbXML schema; making use of IFC entity 'IfcRelSpaceBoundary'
    # This new element is added as child to the earlier created 'gbXML' element
    ifc_global_ids = []

    for ifc_rel_space_boundary in ifc_file.by_type("IfcRelSpaceBoundary"):
        # Make sure a 'SpaceBoundary' is representing an actual element
        ifc_building_element = ifc_rel_space_boundary.RelatedBuildingElement
        if ifc_building_element is None:
            continue

        if ifc_building_element.is_a() in [
            "IfcWall",
            "IfcSlab",
            "IfcRoof",
            "IfcCovering",
        ]:

            # Refer to the relating 'IfcRelAssociatesMaterial' GUID by iterating through IFC entities
            ifc_global_id = ifc_building_element.HasAssociations[0].GlobalId

            # Make use of a list to make sure no same 'Construction' elements are added twice
            if ifc_global_id not in ifc_global_ids:
                ifc_global_ids.append(ifc_global_id)

                # FIXME creates a Construction for every Element. should use Types
                construction = root.createElement("Construction")
                construction.setAttribute("id", fix_xml_cons(ifc_global_id))
                dict_id[fix_xml_cons(ifc_global_id)] = construction

                # Specify analytical properties of the 'Construction' element by iterating through IFC entities

                pset_u_value = (
                    get_pset(
                        # IFC4
                        ifc_building_element,
                        "Pset_WallCommon",
                        prop="ThermalTransmittance",
                    )
                    or get_pset(
                        # IFC4
                        ifc_building_element,
                        "Pset_SlabCommon",
                        prop="ThermalTransmittance",
                    )
                    or get_pset(
                        # IFC4
                        ifc_building_element,
                        "Pset_RoofCommon",
                        prop="ThermalTransmittance",
                    )
                    or get_pset(
                        # IFC4
                        ifc_building_element,
                        "Pset_CoveringCommon",
                        prop="ThermalTransmittance",
                    )
                    or get_pset(
                        # IFC2X3
                        ifc_building_element,
                        "Analytical Properties(Type)",
                        prop="Heat Transfer Coefficient (U)",
                    )
                )
                if pset_u_value:
                    # Building Element could have an overall u-value property rather than layers with thicknesses/r-values
                    u_value = root.createElement("U-value")
                    # FIXME SI units
                    u_value.setAttribute("unit", "WPerSquareMeterK")
                    u_value.appendChild(root.createTextNode(str(pset_u_value)))
                    construction.appendChild(u_value)

                pset_absorptance = get_pset(
                    # IFC2X3, FIXME IFC4 equivalent?
                    ifc_building_element,
                    "Analytical Properties(Type)",
                    prop="Absorptance",
                )
                if pset_absorptance:
                    absorptance = root.createElement("Absorptance")
                    absorptance.setAttribute("unit", "Fraction")
                    absorptance.setAttribute("type", "ExtIR")
                    absorptance.appendChild(root.createTextNode(str(pset_absorptance)))
                    construction.appendChild(absorptance)

                for association in ifc_building_element.HasAssociations:
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

    # NOTE the 'Layer' element of the gbXML schema is not a layer, it is a
    # collection of layers, ie. a layer set
    for ifc_building_element in ifc_file.by_type("IfcBuildingElement"):

        if ifc_building_element.is_a() in [
            "IfcWall",
            "IfcSlab",
            "IfcRoof",
            "IfcCovering",
        ]:

            # Try and catch an Element that is just an Aggregate
            if ifc_building_element.IsDecomposedBy:
                continue

            for association in ifc_building_element.HasAssociations:
                if not association.is_a("IfcRelAssociatesMaterial"):
                    continue
                # FIXME Assumes the IFC element has a Usage
                # FIXME there is a Usage for every instance, so we get a
                # 'Layer' element for every wall, slab and roof, but we don't
                # use Usage attributes. Better to use the layer set for the
                # type (or instance)
                if not association.RelatingMaterial.is_a("IfcMaterialLayerSetUsage"):
                    continue

                # FIXME creates a Layer for every Element. should use Types
                layer = root.createElement("Layer")
                layer.setAttribute("id", fix_xml_layer(association.GlobalId))
                dict_id[fix_xml_layer(association.GlobalId)] = layer

                # Specify the 'IfcMaterialLayer' entity and iterate to each 'IfcMaterial' entity
                for (
                    ifc_material_layer
                ) in association.RelatingMaterial.ForLayerSet.MaterialLayers:
                    material_id = root.createElement("MaterialId")
                    material_id.setAttribute(
                        "materialIdRef", "mat_%d" % ifc_material_layer.id()
                    )
                    layer.appendChild(material_id)

                    dict_id["mat_%d" % ifc_material_layer.id()] = layer

                gbxml.appendChild(layer)

    # NOTE 'Material' element of the gbXML schema is *not* a material, it is a
    # material & thickness combo, ie. a material layer
    ifc_material_layer_ids = []

    for ifc_building_element in ifc_file.by_type("IfcBuildingElement"):

        if ifc_building_element.is_a() in [
            "IfcWall",
            "IfcSlab",
            "IfcRoof",
            "IfcCovering",
        ]:

            # Try and catch an Element that is just an Aggregate
            if ifc_building_element.IsDecomposedBy:
                continue

            for association in ifc_building_element.HasAssociations:
                if not association.is_a("IfcRelAssociatesMaterial"):
                    continue
                # FIXME Assumes the IFC element has a Usage
                if not association.RelatingMaterial.is_a("IfcMaterialLayerSetUsage"):
                    continue

                for (
                    ifc_material_layer
                ) in association.RelatingMaterial.ForLayerSet.MaterialLayers:

                    ifc_material = ifc_material_layer.Material
                    ifc_material_layer_id = ifc_material_layer.id()
                    # Make use of a list to make sure no same 'Materials' elements are added twice
                    if ifc_material_layer_id not in ifc_material_layer_ids:
                        ifc_material_layer_ids.append(ifc_material_layer_id)

                        material = root.createElement("Material")
                        material.setAttribute("id", "mat_%d" % ifc_material_layer_id)
                        dict_id["mat_%d" % ifc_material_layer_id] = material

                        name = root.createElement("Name")
                        name.appendChild(root.createTextNode(ifc_material.Name))
                        material.appendChild(name)

                        thickness = root.createElement("Thickness")
                        # FIXME SI units
                        thickness.setAttribute("unit", "Meters")
                        layer_thickness = ifc_material_layer.LayerThickness
                        thickness.appendChild(
                            root.createTextNode((str(layer_thickness)))
                        )
                        material.appendChild(thickness)

                        r_value = root.createElement("R-value")
                        # FIXME SI units
                        r_value.setAttribute("unit", "SquareMeterKPerW")
                        # NOTE silently sets a default r-value
                        r_value.appendChild(root.createTextNode("0.01"))
                        material.appendChild(r_value)

                        pset_r_value = get_pset(
                            # IFC4
                            ifc_material,
                            "Pset_MaterialEnergy",
                            prop="ThermalConductivityTemperatureDerivative",
                        ) or get_pset(
                            # IFC2X3
                            ifc_material,
                            "Analytical Properties(Type)",
                            prop="Thermal Resistance (R)",
                        )
                        pset_u_value = get_pset(
                            # IFC2X3
                            ifc_material,
                            "Analytical Properties(Type)",
                            prop="Heat Transfer Coefficient (U)",
                        )
                        if pset_u_value and not pset_r_value:
                            pset_r_value = layer_thickness / pset_u_value
                        if pset_r_value:
                            r_value.firstChild.data = str(pset_r_value)

                        gbxml.appendChild(material)

    gbxml.appendChild(create_DocumentHistory(ifc_file, root))
    return root


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

    if not ifc_file.by_type("IfcApplication"):
        program_info = root.createElement("ProgramInfo")
        program_info.setAttribute("id", "IFC_gbXML_Convert")
        document_history.appendChild(program_info)

    for ifc_person in ifc_file.by_type("IfcPerson"):

        created_by = root.createElement("CreatedBy")
        created_by.setAttribute(
            "personId", str(ifc_person.GivenName) + " " + str(ifc_person.FamilyName)
        )
        for ifc_application in ifc_file.by_type("IfcApplication"):
            created_by.setAttribute("programId", ifc_application.ApplicationIdentifier)
        if not ifc_file.by_type("IfcApplication"):
            created_by.setAttribute("programId", "IFC_gbXML_Convert")
        today = datetime.date.today()
        created_by.setAttribute(
            "date", today.strftime("%Y-%m-%dT") + time.strftime("%H:%M:%S")
        )

        document_history.appendChild(created_by)

        person_info = root.createElement("PersonInfo")
        person_info.setAttribute(
            "id", str(ifc_person.GivenName) + " " + str(ifc_person.FamilyName)
        )

        document_history.appendChild(person_info)

    return document_history


# import blenderbim.tool
# ifc_file = blenderbim.tool.Ifc.get()
# root = create_gbxml(ifc_file)
# root.writexml(open("temp.xml", "w"), indent="  ", addindent="  ", newl="\n")

if __name__ == "__main__":
    if not len(sys.argv) == 3:
        print("Usage: " + sys.argv[0] + " input.ifc output.xml")
    else:
        path_ifc = sys.argv[1]
        path_xml = sys.argv[2]
        ifc_file = ifcopenshell.open(path_ifc)

        # Create a new XML file and write all created elements to it
        root = create_gbxml(ifc_file)
        root.writexml(open(path_xml, "w"), indent="  ", addindent="  ", newl="\n")
