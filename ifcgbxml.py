#!/usr/bin/python3

import ifcopenshell.util.placement
import ifcopenshell.util.unit
import numpy as np
import datetime
import time
import sys
import os
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


def get_poly_loop(root, vertices, linear_unit_scale=1.0):
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
            c *= linear_unit_scale
            coordinate = root.createElement("Coordinate")
            coordinate.appendChild(root.createTextNode(str(c)))
            cartesian_point.appendChild(coordinate)
        poly_loop.appendChild(cartesian_point)
    return poly_loop


# valid GUID chars: ^[0123][0-9a-zA-Z_$]{21}$
# valid XML id chars: ^[a-zA-Z_][a-zA-Z0-9_.-]$


def remove_unnecessary_characters(element):
    char_to_replace = {"$": "-", ":": "", " ": "_", "(": "", ")": ""}
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


def fix_xml_bnd(element):
    return "boundary_" + element.replace("$", "-")


def fix_xml_id(element):
    return "id_" + element.replace("$", "-")


# Only fix_xml_name() is used for non-guid data
def fix_xml_name(element):
    return "object_" + remove_unnecessary_characters(element)


def fix_xml_cons(element):
    return "construction_" + element.replace("$", "-")


def fix_xml_layer(element):
    return "layer_" + element.replace("$", "-")


def is_a_boundary_element(ifc_building_element):
    for ifc_class in [
        "IfcWall",
        "IfcColumn",
        "IfcMember",
        "IfcVirtualElement",
        "IfcPlate",
        "IfcCurtainWall",
        "IfcWindow",
        "IfcSlab",
        "IfcRoof",
        "IfcCovering",
    ]:
        if ifc_building_element.is_a(ifc_class):
            return True
    return False


def is_external(ifc_building_element):
    psets = ifcopenshell.util.element.get_psets(
        ifc_building_element, psets_only=True, should_inherit=True
    )
    for pset in psets:
        if "IsExternal" in psets[pset]:
            return True
    return False


def get_thermal_transmittance(ifc_building_element):
    psets = ifcopenshell.util.element.get_psets(
        ifc_building_element, psets_only=True, should_inherit=True
    )
    for pset in psets:
        for prop in psets[pset]:
            if prop == "ThermalTransmittance":
                return psets[pset][prop]
    return get_pset(
        # IFC2X3
        ifc_building_element,
        "Analytical Properties(Type)",
        prop="Heat Transfer Coefficient (U)",
    )


def get_parent_boundary(ifc_rel_space_boundary):
    ifc_building_element = ifc_rel_space_boundary.RelatedBuildingElement
    if hasattr(ifc_rel_space_boundary, "ParentBoundary"):
        # IFC4
        return ifc_rel_space_boundary.ParentBoundary
    elif (
        ifc_building_element.FillsVoids
        and ifc_building_element.FillsVoids.RelatingOpeningElement.VoidsElements
    ):
        # IFC2X3
        for ifc_boundary in (
            ifc_building_element.FillsVoids[0]
            .RelatingOpeningElement.VoidsElements[0]
            .RelatingBuildingElement.ProvidesBoundaries
        ):
            if ifc_boundary.RelatingSpace == ifc_rel_space_boundary.RelatingSpace:
                return ifc_boundary


def get_material_layer_set(ifc_building_element):
    for association in ifc_building_element.HasAssociations:
        if association.is_a("IfcRelAssociatesMaterial"):
            if association.RelatingMaterial.is_a("IfcMaterialLayerSet"):
                return association.RelatingMaterial
            elif association.RelatingMaterial.is_a("IfcMaterialLayerSetUsage"):
                return association.RelatingMaterial.ForLayerSet


def create_gbxml(ifc_file):
    """Process an IfcOpenShell file object and return a minidom document in gbXML format"""
    root = minidom.Document()

    gbxml = root.createElement("gbXML")
    root.appendChild(gbxml)

    gbxml.setAttribute("xmlns", "http://www.gbxml.org/schema")
    gbxml.setAttribute("temperatureUnit", "C")
    gbxml.setAttribute("lengthUnit", "Meters")
    gbxml.setAttribute("areaUnit", "SquareMeters")
    gbxml.setAttribute("volumeUnit", "CubicMeters")
    gbxml.setAttribute("useSIUnitsForResults", "true")
    gbxml.setAttribute("version", "0.37")

    linear_unit_scale = ifcopenshell.util.unit.calculate_unit_scale(ifc_file)

    area_unit_scale = 1.0
    area_unit = ifcopenshell.util.unit.get_project_unit(ifc_file, "AREAUNIT")
    while area_unit.is_a("IfcConversionBasedUnit"):
        area_unit_scale *= area_unit.ConversionFactor.ValueComponent.wrappedValue
        area_unit = area_unit.ConversionFactor.UnitComponent
    if area_unit.is_a("IfcSIUnit"):
        area_unit_scale *= ifcopenshell.util.unit.get_prefix_multiplier(
            area_unit.Prefix
        )

    volume_unit_scale = 1.0
    volume_unit = ifcopenshell.util.unit.get_project_unit(ifc_file, "VOLUMEUNIT")
    while volume_unit.is_a("IfcConversionBasedUnit"):
        volume_unit_scale *= volume_unit.ConversionFactor.ValueComponent.wrappedValue
        volume_unit = volume_unit.ConversionFactor.UnitComponent
    if volume_unit.is_a("IfcSIUnit"):
        volume_unit_scale *= ifcopenshell.util.unit.get_prefix_multiplier(
            volume_unit.Prefix
        )

    # NOTE crude check for imperial units
    imperial_units = False
    if ifcopenshell.util.unit.get_project_unit(ifc_file, "LENGTHUNIT").Name in [
        "inch",
        "foot",
        "yard",
    ]:
        imperial_units = True

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
                name.appendChild(
                    root.createTextNode(ifc_building_storey.Name or "Unnamed")
                )
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

                    space.setAttribute(
                        "buildingStoreyIdRef",
                        fix_xml_stry(ifc_space.Decomposes[0].RelatingObject.GlobalId),
                    )

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
                        area = root.createElement("Area")
                        area.setAttribute("unit", "SquareMeters")
                        area.appendChild(
                            root.createTextNode(str(pset_area * area_unit_scale))
                        )
                        space.appendChild(area)

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
                        volume = root.createElement("Volume")
                        volume.setAttribute("unit", "CubicMeters")
                        volume.appendChild(
                            root.createTextNode(str(pset_volume * volume_unit_scale))
                        )
                        space.appendChild(volume)

                    name = root.createElement("Name")
                    name.appendChild(root.createTextNode(ifc_space.Name or "Unnamed"))
                    space.appendChild(name)

                    # Specify the 'SpaceBoundary' element of the gbXML schema; making use of IFC entity 'IfcSpace'
                    # This new element is added as child to the earlier created 'Space' element
                    for ifc_rel_space_boundary in ifc_space.BoundedBy:

                        if get_parent_boundary(ifc_rel_space_boundary):
                            continue

                        # Make sure a 'SpaceBoundary' is representing an actual element
                        if ifc_rel_space_boundary.RelatedBuildingElement == None:
                            continue
                        if (
                            ifc_rel_space_boundary.ConnectionGeometry.SurfaceOnRelatingElement
                            == None
                        ):
                            continue
                        # require 2nd Level boundaries
                        if not (
                            ifc_rel_space_boundary.Name == "2ndLevel"
                            or ifc_rel_space_boundary.is_a(
                                "IfcRelSpaceBoundary2ndLevel"
                            )
                        ):
                            continue

                        vertices = get_boundary_vertices(ifc_rel_space_boundary)

                        ifc_building_element = (
                            ifc_rel_space_boundary.RelatedBuildingElement
                        )

                        # Create 'SpaceBoundary' elements for the following building elements
                        if is_a_boundary_element(ifc_building_element):

                            space_boundary = root.createElement("SpaceBoundary")
                            space_boundary.setAttribute("isSecondLevelBoundary", "true")

                            space_boundary.setAttribute(
                                "surfaceIdRef",
                                fix_xml_id(ifc_rel_space_boundary.GlobalId),
                            )

                            space.appendChild(space_boundary)

                            planar_geometry = root.createElement("PlanarGeometry")
                            space_boundary.appendChild(planar_geometry)

                            planar_geometry.appendChild(
                                get_poly_loop(
                                    root, vertices, linear_unit_scale=linear_unit_scale
                                )
                            )

    # Specify the 'Surface' element of the gbXML schema; making use of IFC entity 'IfcRelSpaceBoundary'
    # This new element is added as child to the earlier created 'Campus' element
    for ifc_rel_space_boundary in ifc_file.by_type("IfcRelSpaceBoundary"):

        ifc_building_element = ifc_rel_space_boundary.RelatedBuildingElement

        # Make sure a 'SpaceBoundary' is representing an actual element
        if ifc_building_element == None:
            continue
        if ifc_rel_space_boundary.ConnectionGeometry.SurfaceOnRelatingElement == None:
            continue
        if ifc_rel_space_boundary.RelatingSpace == None:
            continue
        # require 2nd Level boundaries
        if not (
            ifc_rel_space_boundary.Name == "2ndLevel"
            or ifc_rel_space_boundary.is_a("IfcRelSpaceBoundary2ndLevel")
        ):
            continue
        if get_parent_boundary(ifc_rel_space_boundary):
            continue

        vertices = get_boundary_vertices(ifc_rel_space_boundary)

        # Specify each 'Surface' element and set 'SurfaceType' attributes
        if is_a_boundary_element(ifc_building_element):

            surface = root.createElement("Surface")
            surface.setAttribute("id", fix_xml_id(ifc_rel_space_boundary.GlobalId))
            dict_id[fix_xml_id(ifc_rel_space_boundary.GlobalId)] = surface

            if (
                hasattr(ifc_building_element, "IsTypedBy")
                and ifc_building_element.IsTypedBy
            ):
                ifc_building_element = ifc_building_element.IsTypedBy[0].RelatingType

            # Valid surfaceType:
            # InteriorWall, ExteriorWall, Roof, InteriorFloor, ExposedFloor,
            # Shade, UndergroundWall, UndergroundSlab, Ceiling, Air,
            # UndergroundCeiling, RaisedFloor, SlabOnGrade, FreestandingColumn,
            # EmbeddedColumn, Undefined

            if ifc_building_element.is_a("IfcCovering") or ifc_building_element.is_a(
                "IfcCoveringType"
            ):
                # NOTE assumes all coverings with a related space boundary are ceilings
                surface.setAttribute("surfaceType", "Ceiling")

            elif ifc_building_element.is_a("IfcRoof") or ifc_building_element.is_a(
                "IfcRoofType"
            ):
                surface.setAttribute("surfaceType", "Roof")
                surface.setAttribute("exposedToSun", "true")

            elif ifc_building_element.is_a("IfcColumn") or ifc_building_element.is_a(
                "IfcColumnType"
            ):
                surface.setAttribute("surfaceType", "EmbeddedColumn")

            elif ifc_building_element.is_a("IfcSlab") or ifc_building_element.is_a(
                "IfcSlabType"
            ):
                if (
                    ifc_rel_space_boundary.InternalOrExternalBoundary
                    == "EXTERNAL_EARTH"
                ) or get_pset(
                    ifc_building_element, "Pset_SlabCommon", prop="IsExternal"
                ):
                    surface.setAttribute("surfaceType", "SlabOnGrade")
                else:
                    surface.setAttribute("surfaceType", "InteriorFloor")

            elif ifc_building_element.is_a("IfcWall") or ifc_building_element.is_a(
                "IfcWallType"
            ):
                if (
                    ifc_rel_space_boundary.InternalOrExternalBoundary == "EXTERNAL"
                    or get_pset(
                        ifc_building_element, "Pset_WallCommon", prop="IsExternal"
                    )
                ):
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

            elif ifc_building_element.is_a(
                "IfcCurtainWall"
            ) or ifc_building_element.is_a("IfcCurtainWallType"):
                surface.setAttribute("surfaceType", "ExteriorWall")
                surface.setAttribute("exposedToSun", "true")

            elif ifc_building_element.is_a("IfcVirtualElement"):
                surface.setAttribute("surfaceType", "Air")

            else:
                surface.setAttribute("surfaceType", "Undefined")
                if is_external(ifc_building_element):
                    surface.setAttribute("exposedToSun", "true")

            ifc_material_layer_set = get_material_layer_set(ifc_building_element)

            if ifc_material_layer_set:
                surface.setAttribute(
                    "constructionIdRef",
                    fix_xml_cons(str(ifc_material_layer_set.id())),
                )
            else:
                # elements without a layer set may have u-value property
                surface.setAttribute(
                    "constructionIdRef",
                    fix_xml_cons(str(ifc_building_element.id())),
                )

            name = root.createElement("Name")
            name.appendChild(
                root.createTextNode(fix_xml_bnd(ifc_rel_space_boundary.GlobalId))
            )

            surface.appendChild(name)

            adjacent_space_id = root.createElement("AdjacentSpaceId")

            adjacent_space_id.setAttribute(
                "spaceIdRef", fix_xml_spc(ifc_rel_space_boundary.RelatingSpace.GlobalId)
            )
            surface.appendChild(adjacent_space_id)

            planar_geometry = root.createElement("PlanarGeometry")
            surface.appendChild(planar_geometry)

            planar_geometry.appendChild(
                get_poly_loop(root, vertices, linear_unit_scale=linear_unit_scale)
            )

            cad_object_id = root.createElement("CADObjectId")
            cad_object_id.appendChild(
                root.createTextNode(fix_xml_bnd(ifc_rel_space_boundary.GlobalId))
            )
            surface.appendChild(cad_object_id)

            campus.appendChild(surface)

    # Specify the 'Opening' element of the gbXML schema; Windows and Doors
    # This new element is added as child to an earlier created 'Surface' element
    opening_id = 1
    for ifc_rel_space_boundary in ifc_file.by_type("IfcRelSpaceBoundary"):

        ifc_building_element = ifc_rel_space_boundary.RelatedBuildingElement
        # Make sure a 'SpaceBoundary' is representing an actual element
        if ifc_building_element == None:
            continue
        if ifc_rel_space_boundary.ConnectionGeometry.SurfaceOnRelatingElement == None:
            continue

        vertices = get_boundary_vertices(ifc_rel_space_boundary)

        if ifc_building_element.is_a() in [
            "IfcWindow",
            "IfcDoor",
        ]:

            ifc_parent_boundary = get_parent_boundary(ifc_rel_space_boundary)

            if not ifc_parent_boundary:
                continue

            if (
                hasattr(ifc_building_element, "IsTypedBy")
                and ifc_building_element.IsTypedBy
            ):
                ifc_building_element = ifc_building_element.IsTypedBy[0].RelatingType

            opening = root.createElement("Opening")
            dict_id[fix_xml_id(ifc_rel_space_boundary.GlobalId)] = opening

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

            planar_geometry.appendChild(
                get_poly_loop(root, vertices, linear_unit_scale=linear_unit_scale)
            )

            name = root.createElement("Name")
            name.appendChild(
                root.createTextNode(
                    fix_xml_name(ifc_building_element.Name or "Unnamed")
                )
            )
            opening.appendChild(name)

            cad_object_id = root.createElement("CADObjectId")
            cad_object_id.appendChild(
                root.createTextNode(
                    fix_xml_name(ifc_building_element.Name or "Unnamed")
                )
            )
            opening.appendChild(cad_object_id)

            if fix_xml_id(ifc_parent_boundary.GlobalId) in dict_id:
                surface = dict_id[fix_xml_id(ifc_parent_boundary.GlobalId)]
                surface.appendChild(opening)

    ifc_building_element_guids = []
    # Specify the 'WindowType' element of the gbXML schema; making use of IFC entity 'IfcWindow'
    # This new element is added as child to the earlier created 'gbXML' element
    for ifc_building_element in ifc_file.by_type("IfcWindow") + ifc_file.by_type(
        "IfcDoor"
    ):

        if (
            hasattr(ifc_building_element, "IsTypedBy")
            and ifc_building_element.IsTypedBy
        ):
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
                root.createTextNode(
                    fix_xml_name(ifc_building_element.Name or "Unnamed")
                )
            )
            window_type.appendChild(name)

            description = root.createElement("Description")
            description.appendChild(
                root.createTextNode(
                    fix_xml_name(ifc_building_element.Name or "Unnamed")
                )
            )
            window_type.appendChild(description)

            u_value = root.createElement("U-value")
            u_value.setAttribute("unit", "WPerSquareMeterK")
            u_value.appendChild(root.createTextNode("10.0"))

            pset_u_value = get_thermal_transmittance(ifc_building_element)

            if pset_u_value:
                if imperial_units:
                    pset_u_value *= 5.678
                u_value.firstChild.data = str(pset_u_value)
                window_type.appendChild(u_value)

            solar_heat_gain_coeff = root.createElement("SolarHeatGainCoeff")
            solar_heat_gain_coeff.setAttribute("unit", "Fraction")
            solar_heat_gain_coeff.appendChild(root.createTextNode("1.0"))

            # FIXME set default, IFC4?
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
    ifc_ids = []
    ifc_material_layer_ids = []

    for ifc_rel_space_boundary in ifc_file.by_type("IfcRelSpaceBoundary"):
        # Make sure a 'SpaceBoundary' is representing an actual element
        ifc_building_element = ifc_rel_space_boundary.RelatedBuildingElement
        if ifc_building_element is None:
            continue
        if get_parent_boundary(ifc_rel_space_boundary):
            continue

        if is_a_boundary_element(ifc_building_element):
            if (
                hasattr(ifc_building_element, "IsTypedBy")
                and ifc_building_element.IsTypedBy
            ):
                ifc_building_element = ifc_building_element.IsTypedBy[0].RelatingType

            ifc_material_layer_set = get_material_layer_set(ifc_building_element)

            if ifc_material_layer_set:
                ifc_id = str(ifc_material_layer_set.id())
                construction_name = ifc_material_layer_set.LayerSetName or "Unnamed"
            else:
                # elements without a layer set may have u-value property
                ifc_id = str(ifc_building_element.id())
                construction_name = ifc_building_element.Name or "Unnamed"

            # Make use of a list to make sure no same 'Construction' elements are added twice
            if ifc_id not in ifc_ids:
                ifc_ids.append(ifc_id)

                construction = root.createElement("Construction")
                construction.setAttribute("id", fix_xml_cons(ifc_id))
                dict_id[fix_xml_cons(ifc_id)] = construction

                pset_u_value = get_thermal_transmittance(ifc_building_element)

                if pset_u_value:
                    # Building Element could have an overall u-value property rather than layers with thicknesses/r-values
                    u_value = root.createElement("U-value")
                    if imperial_units:
                        pset_u_value *= 5.678
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

                name = root.createElement("Name")
                name.appendChild(root.createTextNode(construction_name))
                construction.appendChild(name)

                gbxml.appendChild(construction)

                if not ifc_material_layer_set:
                    # elements without a layer set may have u-value property
                    continue

                layer_id = root.createElement("LayerId")
                layer_id.setAttribute("layerIdRef", fix_xml_layer(ifc_id))
                construction.appendChild(layer_id)

                # NOTE the 'Layer' element of the gbXML schema is not a layer, it is a
                # collection of layers, ie. a layer set
                layer = root.createElement("Layer")
                layer.setAttribute("id", fix_xml_layer(ifc_id))
                dict_id[fix_xml_layer(ifc_id)] = layer

                for ifc_material_layer in ifc_material_layer_set.MaterialLayers:
                    material_id = root.createElement("MaterialId")
                    material_id.setAttribute(
                        "materialIdRef", "mat_%d" % ifc_material_layer.id()
                    )
                    layer.appendChild(material_id)

                    dict_id["mat_%d" % ifc_material_layer.id()] = layer

                gbxml.appendChild(layer)

                for ifc_material_layer in ifc_material_layer_set.MaterialLayers:

                    ifc_material = ifc_material_layer.Material
                    ifc_material_layer_id = ifc_material_layer.id()
                    # Make use of a list to make sure no same 'Materials' elements are added twice
                    if ifc_material_layer_id not in ifc_material_layer_ids:
                        ifc_material_layer_ids.append(ifc_material_layer_id)

                        material = root.createElement("Material")
                        material.setAttribute("id", "mat_%d" % ifc_material_layer_id)
                        dict_id["mat_%d" % ifc_material_layer_id] = material

                        name = root.createElement("Name")
                        name.appendChild(
                            root.createTextNode(ifc_material.Name or "Unnamed")
                        )
                        material.appendChild(name)

                        thickness = root.createElement("Thickness")
                        thickness.setAttribute("unit", "Meters")
                        layer_thickness = (
                            ifc_material_layer.LayerThickness * linear_unit_scale
                        )
                        thickness.appendChild(
                            root.createTextNode((str(layer_thickness)))
                        )
                        material.appendChild(thickness)

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
                        if pset_r_value:
                            if imperial_units:
                                pset_r_value /= 5.678
                            r_value = root.createElement("R-value")
                            r_value.setAttribute("unit", "SquareMeterKPerW")
                            r_value.appendChild(root.createTextNode(str(pset_r_value)))
                            material.appendChild(r_value)
                        elif pset_u_value:
                            r_value = root.createElement("R-value")
                            r_value.setAttribute("unit", "SquareMeterKPerW")
                            r_value.appendChild(
                                root.createTextNode(str(layer_thickness / pset_u_value))
                            )
                            material.appendChild(r_value)

                        gbxml.appendChild(material)

    gbxml.appendChild(create_DocumentHistory(ifc_file, root))
    return root


def create_DocumentHistory(ifc_file, root):
    """Specify the 'DocumentHistory' element of the gbXML schema; making use of IFC entity 'IfcApplication' and 'IfcPerson'"""
    document_history = root.createElement("DocumentHistory")

    program_info = root.createElement("ProgramInfo")
    program_info.setAttribute("id", "IFC_gbXML_Convert")
    document_history.appendChild(program_info)

    person_info = root.createElement("PersonInfo")
    person_info.setAttribute("id", remove_unnecessary_characters(os.getlogin()))
    document_history.appendChild(person_info)

    created_by = root.createElement("CreatedBy")
    created_by.setAttribute("personId", remove_unnecessary_characters(os.getlogin()))
    created_by.setAttribute("programId", "IFC_gbXML_Convert")
    today = datetime.date.today()
    created_by.setAttribute(
        "date", today.strftime("%Y-%m-%dT") + time.strftime("%H:%M:%S")
    )
    document_history.appendChild(created_by)

    return document_history


# import blenderbim.tool
# root = create_gbxml(blenderbim.tool.Ifc.get())
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
