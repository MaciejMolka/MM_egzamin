# -*- coding: utf-8 -*-

import arcpy
from math import sqrt, acos, pi
import numpy as np
import pandas as pd

def tworzenie_warstwy_punktowej(pnts):
    arcpy.CreateFeatureclass_management("C:/Studia/PPG_II/Egzamin/Wyniki", "Concave_Hull.shp", "POINT")
    cur = arcpy.da.InsertCursor("C:/Studia/PPG_II/Egzamin/Wyniki/Concave_Hull.shp", ["SHAPE@"])

    punkty = []
    for i in pnts:
        for j in i:
            punkty.append(j)

    for row in punkty:
        cur.insertRow([row])

def concave_hull(warstwa_wejsciowa, warstwa_wynikowa, k_pocz=3, pole="", puste=True):
    try:
        import arcpy
        import itertools
        import math
        import os
        import sys
        import traceback
        import string

        arcpy.overwriteOutput = True

        # Functions that consolidate reuable actions
        #

        # Function to return an OID list for k nearest eligible neighbours of a feature
        def kNeighbours(k, oid, pDict, excludeList=[]):
            hypotList = [math.hypot(pDict[oid][0] - pDict[id][0], pDict[oid][1] - pDict[id][1]) for id in pDict.keys()
                         if
                         id <> oid and id not in excludeList]
            hypotList.sort()
            hypotList = hypotList[0:k]
            oidList = [id for id in pDict.keys() if math.hypot(pDict[oid][0] - pDict[id][0], pDict[oid][1] - pDict[id][
                1]) in hypotList and id <> oid and id not in excludeList]
            return oidList

        # Function to rotate a point about another point, returning a list [X,Y]
        def RotateXY(x, y, xc=0, yc=0, angle=0):
            x = x - xc
            y = y - yc
            xr = (x * math.cos(angle)) - (y * math.sin(angle)) + xc
            yr = (x * math.sin(angle)) + (y * math.cos(angle)) + yc
            return [xr, yr]

        # Function finding the feature OID at the rightmost angle from an origin OID, with respect to an input angle
        def Rightmost(oid, angle, pDict, oidList):
            origxyList = [pDict[id] for id in pDict.keys() if id in oidList]
            rotxyList = []
            for p in range(len(origxyList)):
                rotxyList.append(RotateXY(origxyList[p][0], origxyList[p][1], pDict[oid][0], pDict[oid][1], angle))
            minATAN = min([math.atan2((xy[1] - pDict[oid][1]), (xy[0] - pDict[oid][0])) for xy in rotxyList])
            rightmostIndex = rotxyList.index(
                [xy for xy in rotxyList if math.atan2((xy[1] - pDict[oid][1]), (xy[0] - pDict[oid][0])) == minATAN][0])
            return oidList[rightmostIndex]

        # Function to detect single-part polyline self-intersection
        def selfIntersects(polyline):
            lList = []
            selfIntersects = False
            for n in range(0, len(line.getPart(0)) - 1):
                lList.append(arcpy.Polyline(arcpy.Array([line.getPart(0)[n], line.getPart(0)[n + 1]])))
            for pair in itertools.product(lList, repeat=2):
                if pair[0].crosses(pair[1]):
                    selfIntersects = True
                    break
            return selfIntersects

        # Function to construct the Hull
        def createHull(pDict, outCaseField, lastValue, kStart, dictCount, includeNull):
            # Value of k must result in enclosing all data points; create condition flag
            enclosesPoints = False
            notNullGeometry = False
            k = kStart

            if dictCount > 1:
                pList = [arcpy.Point(xy[0], xy[1]) for xy in pDict.values()]
                mPoint = arcpy.Multipoint(arcpy.Array(pList), sR)
                minY = min([xy[1] for xy in pDict.values()])

                while not enclosesPoints and k <= 30:
                    arcpy.AddMessage("Finding hull for k = " + str(k))
                    # Find start point (lowest Y value)
                    startOID = [id for id in pDict.keys() if pDict[id][1] == minY][0]
                    # Select the next point (rightmost turn from horizontal, from start point)
                    kOIDList = kNeighbours(k, startOID, pDict, [])
                    minATAN = min(
                        [math.atan2(pDict[id][1] - pDict[startOID][1], pDict[id][0] - pDict[startOID][0]) for id in
                         kOIDList])
                    nextOID = [id for id in kOIDList if
                               math.atan2(pDict[id][1] - pDict[startOID][1],
                                          pDict[id][0] - pDict[startOID][0]) == minATAN][
                        0]
                    # Initialise the boundary array
                    bArray = arcpy.Array(arcpy.Point(pDict[startOID][0], pDict[startOID][1]))
                    bArray.add(arcpy.Point(pDict[nextOID][0], pDict[nextOID][1]))
                    # Initialise current segment lists
                    currentOID = nextOID
                    prevOID = startOID
                    # Initialise list to be excluded from candidate consideration (start point handled additionally later)
                    excludeList = [startOID, nextOID]
                    # Build the boundary array - taking the closest rightmost point that does not cause a self-intersection.
                    steps = 2
                    while currentOID <> startOID and len(pDict) <> len(excludeList):
                        try:
                            angle = math.atan2((pDict[currentOID][1] - pDict[prevOID][1]),
                                               (pDict[currentOID][0] - pDict[prevOID][0]))
                            oidList = kNeighbours(k, currentOID, pDict, excludeList)
                            nextOID = Rightmost(currentOID, 0 - angle, pDict, oidList)
                            pcArray = arcpy.Array([arcpy.Point(pDict[currentOID][0], pDict[currentOID][1]), \
                                                   arcpy.Point(pDict[nextOID][0], pDict[nextOID][1])])
                            while arcpy.Polyline(bArray, sR).crosses(arcpy.Polyline(pcArray, sR)) and len(oidList) > 0:
                                # arcpy.AddMessage("Rightmost point from " + str(currentOID) + " : " + str(nextOID) + " causes self intersection - selecting again")
                                excludeList.append(nextOID)
                                oidList.remove(nextOID)
                                oidList = kNeighbours(k, currentOID, pDict, excludeList)
                                if len(oidList) > 0:
                                    nextOID = Rightmost(currentOID, 0 - angle, pDict, oidList)
                                    # arcpy.AddMessage("nextOID candidate: " + str(nextOID))
                                    pcArray = arcpy.Array([arcpy.Point(pDict[currentOID][0], pDict[currentOID][1]), \
                                                           arcpy.Point(pDict[nextOID][0], pDict[nextOID][1])])
                            bArray.add(arcpy.Point(pDict[nextOID][0], pDict[nextOID][1]))
                            prevOID = currentOID
                            currentOID = nextOID
                            excludeList.append(currentOID)
                            # arcpy.AddMessage("CurrentOID = " + str(currentOID))
                            steps += 1
                            if steps == 4:
                                excludeList.remove(startOID)
                        except ValueError:
                            arcpy.AddMessage(
                                "Zero reachable nearest neighbours at " + str(
                                    pDict[currentOID]) + " , expanding search")
                            break
                    # Close the boundary and test for enclosure
                    bArray.add(arcpy.Point(pDict[startOID][0], pDict[startOID][1]))
                    pPoly = arcpy.Polygon(bArray, sR)
                    if pPoly.length == 0:
                        break
                    else:
                        notNullGeometry = True
                    if mPoint.within(arcpy.Polygon(bArray, sR)):
                        enclosesPoints = True
                    else:
                        arcpy.AddMessage("Hull does not enclose data, incrementing k")
                        k += 1
                #
                if not mPoint.within(arcpy.Polygon(bArray, sR)):
                    arcpy.AddWarning("Hull does not enclose data - probable cause is outlier points")

            # Insert the Polygons
            if (notNullGeometry and includeNull == False) or includeNull:
                if outCaseField > " ":
                    insFields = [outCaseField, "POINT_CNT", "ENCLOSED", "SHAPE@"]
                else:
                    insFields = ["POINT_CNT", "ENCLOSED", "SHAPE@"]
                rows = arcpy.da.InsertCursor(outFC, insFields)
                row = []
                if outCaseField > " ":
                    row.append(lastValue)
                row.append(dictCount)
                if notNullGeometry:
                    row.append(enclosesPoints)
                    row.append(arcpy.Polygon(bArray, sR))
                else:
                    row.append(-1)
                    row.append(None)
                rows.insertRow(row)
                del row
                del rows
            elif outCaseField > " ":
                arcpy.AddMessage("\nExcluded Null Geometry for case value " + str(lastValue) + "!")
            else:
                arcpy.AddMessage("\nExcluded Null Geometry!")

        # Main Body of the program.
        #
        #

        # Get the input feature class or layer
        inPoints = warstwa_wejsciowa
        inDesc = arcpy.Describe(inPoints)
        inPath = os.path.dirname(inDesc.CatalogPath)
        sR = inDesc.spatialReference

        # Get k
        k = int(k_pocz)
        kStart = k

        # Get output Feature Class
        outFC = warstwa_wynikowa
        outPath = os.path.dirname(outFC)
        outName = os.path.basename(outFC)

        # Get case field and ensure it is valid
        caseField = pole
        if caseField > " ":
            fields = inDesc.fields
            for field in fields:
                # Check the case field type
                if field.name == caseField:
                    caseFieldType = field.type
                    if caseFieldType not in ["SmallInteger", "Integer", "Single", "Double", "String", "Date"]:
                        arcpy.AddMessage(
                            "\nThe Case Field named " + caseField + " is not a valid case field type!  The Case Field will be ignored!\n")
                        caseField = " "
                    else:
                        if caseFieldType in ["SmallInteger", "Integer", "Single", "Double"]:
                            caseFieldLength = 0
                            caseFieldScale = field.scale
                            caseFieldPrecision = field.precision
                        elif caseFieldType == "String":
                            caseFieldLength = field.length
                            caseFieldScale = 0
                            caseFieldPrecision = 0
                        else:
                            caseFieldLength = 0
                            caseFieldScale = 0
                            caseFieldPrecision = 0

        # Define an output case field name that is compliant with the output feature class
        outCaseField = str.upper(str(caseField))
        if outCaseField == "ENCLOSED":
            outCaseField = "ENCLOSED1"
        if outCaseField == "POINT_CNT":
            outCaseField = "POINT_CNT1"
        if outFC.split(".")[-1] in ("shp", "dbf"):
            outCaseField = outCaseField[0:10]  # field names in the output are limited to 10 charaters!

        # Get Include Null Geometry Feature flag
        includeNull = puste

        # Some housekeeping
        inDesc = arcpy.Describe(inPoints)
        sR = inDesc.spatialReference
        arcpy.env.OutputCoordinateSystem = sR
        oidName = str(inDesc.OIDFieldName)
        if inDesc.dataType == "FeatureClass":
            inPoints = arcpy.MakeFeatureLayer_management(inPoints)

        # Create the output
        arcpy.AddMessage("\nCreating Feature Class...")
        if '.SHP' in outName.upper():
            outName = outName[:-4]
        arcpy.AddMessage(outPath + "; " + outName)
        outFC = arcpy.CreateFeatureclass_management(outPath, outName, "POLYGON", "#", "#", "#", sR).getOutput(0)
        if caseField > " ":
            if caseFieldType in ["SmallInteger", "Integer", "Single", "Double"]:
                arcpy.AddField_management(outFC, outCaseField, caseFieldType, str(caseFieldScale),
                                          str(caseFieldPrecision))
            elif caseFieldType == "String":
                arcpy.AddField_management(outFC, outCaseField, caseFieldType, "", "", str(caseFieldLength))
            else:
                arcpy.AddField_management(outFC, outCaseField, caseFieldType)
        arcpy.AddField_management(outFC, "POINT_CNT", "Long")
        arcpy.AddField_management(outFC, "ENCLOSED", "SmallInteger")

        # Build required data structures
        arcpy.AddMessage("\nCreating data structures...")
        rowCount = 0
        caseCount = 0
        dictCount = 0
        pDict = {}  # dictionary keyed on oid with [X,Y] list values, no duplicate points
        if caseField > " ":
            fields = [caseField, 'OID@', 'SHAPE@X', 'SHAPE@Y']
            valueDict = {}
            with arcpy.da.SearchCursor(inPoints, fields) as searchRows:
                for searchRow in searchRows:
                    keyValue = searchRow[0]
                    if not keyValue in valueDict:
                        # assign a new keyValue entry to the dictionary storing a list of the first NumberField value and 1 for the first record counter value
                        valueDict[keyValue] = [[searchRow[1], searchRow[2], searchRow[3]]]
                    # Sum the last summary of NumberField value with the current record and increment the record count when keyvalue is already in the dictionary
                    else:
                        valueDict[keyValue].append([searchRow[1], searchRow[2], searchRow[3]])
            for lastValue in sorted(valueDict):
                caseCount += 1
                for p in valueDict[lastValue]:
                    rowCount += 1
                    # Continue processing the current point subset.
                    if [p[1], p[2]] not in pDict.values():
                        pDict[p[0]] = [p[1], p[2]]
                        dictCount += 1
                createHull(pDict, outCaseField, lastValue, kStart, dictCount, includeNull)
                # Reset variables for processing the next point subset.
                pDict = {}
                dictCount = 0
        else:
            fields = ['OID@', 'SHAPE@X', 'SHAPE@Y']
            for p in arcpy.da.SearchCursor(inPoints, fields):
                rowCount += 1
                if [p[1], p[2]] not in pDict.values():
                    pDict[p[0]] = [p[1], p[2]]
                    dictCount += 1
                    lastValue = 0
            # Final create hull call and wrap up of the program's massaging
            createHull(pDict, outCaseField, lastValue, kStart, dictCount, includeNull)
        arcpy.AddMessage("\n" + str(rowCount) + " points processed.  " + str(caseCount) + " case value(s) processed.")
        if caseField == " " and arcpy.GetParameterAsText(3) > " ":
            arcpy.AddMessage(
                "\nThe Case Field named " + arcpy.GetParameterAsText(
                    3) + " was not a valid field type and was ignored!")
        arcpy.AddMessage("\nFinished")


    # Error handling
    except:
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]
        pymsg = "PYTHON ERRORS:\nTraceback Info:\n" + tbinfo + "\nError Info:\n    " + \
                str(sys.exc_type) + ": " + str(sys.exc_value) + "\n"
        arcpy.AddError(pymsg)

        msgs = "GP ERRORS:\n" + arcpy.GetMessages(2) + "\n"
        arcpy.AddError(msgs)

    arcpy.CheckGeometry_management(warstwa_wynikowa, "C:/Studia/PPG_II/Egzamin/Wyniki/spr")

    if arcpy.GetCount_management("C:/Studia/PPG_II/Egzamin/Wyniki/spr")[0] == "2":
        return "Brak"
    else:
        return polygon_do_linii(warstwa_wynikowa)

def geometria_czesci(part, obiekt, sprawdzenie, minimalne_geometrie):
    nr_punktu = -1
    wyniki = []
    for i in part[0:-1]:
        nr_punktu += 1
        punkt_poprzedni = part[0:-1][nr_punktu - 1]
        punkt_srodkowy = part[0:-1][nr_punktu]
        if nr_punktu != len(part[0:-1]) - 1:
            punkt_nastepny = part[0:-1][nr_punktu + 1]
        else:
            punkt_nastepny = part[0:-1][0]

        length_in = segment_length(punkt_poprzedni, punkt_srodkowy)
        length_out = segment_length(punkt_nastepny, punkt_srodkowy)
        angle_in = vertex_angle(punkt_poprzedni, punkt_srodkowy, punkt_nastepny, obiekt, sprawdzenie)

        min_odleglosci = []
        for i in minimalne_geometrie:
            if i == "Brak":
                min_odleglosci.append("Nieokreslone")
            else:
                min_odleglosci.append(deflection(punkt_srodkowy, i))

        wyniki.append(
            [nr_punktu, length_in, length_out, angle_in, min_odleglosci[0], min_odleglosci[1], min_odleglosci[2],
             min_odleglosci[3], min_odleglosci[4], min_odleglosci[5]])

    return wyniki

def geometria_obiektu(parts, obiekt, minimalne_geometrie):
    iterator = 0
    wynik = []
    for i in parts:
        iterator += 1
        if iterator == 1:
            wyniki = geometria_czesci(i, obiekt, 'normalnie', minimalne_geometrie)
            wynik.append(wyniki)
        else:
            anArray = arcpy.Array()
            for j in i:
                anArray.add(j)
            nowy_obiekt = arcpy.Polygon(anArray)
            wyniki = geometria_czesci(i, nowy_obiekt, 'odwrotnosc', minimalne_geometrie)
            wynik.append(wyniki)

    wynik_ostateczny = []
    for i in wynik:
        for j in i:
            wynik_ostateczny.append(j)

    return wynik_ostateczny

def polygon_do_linii(polygon):
    arcpy.PolygonToLine_management(polygon, polygon[0:-4] + "_linia.shp")

    linie = []
    for row in arcpy.da.SearchCursor(polygon[0:-4] + "_linia.shp", ["SHAPE@"]):
        linie.append(row[0])

    return linie[0]

def deflection(wierzcholek, linia):
    return linia.distanceTo(wierzcholek)

def intersect(geometria1, geometria2):
    if not geometria1.disjoint(geometria2):
        return True
    else:
        return False

def is_node(wierzcholek, id_budynku, obiekty):
    expression = "gmlID <> " + "'" + id_budynku + "'"

    ilosc_sasiadow = 0
    for row in arcpy.da.SearchCursor(obiekty, ["OID@", "SHAPE@", "gmlId"], where_clause=expression):
        for part in row[1]:
            polygon = arcpy.Polygon(part)
            if intersect(wierzcholek, polygon):
                ilosc_sasiadow = ilosc_sasiadow + 1

    return ilosc_sasiadow

def minimal_geometry(feature):
    arcpy.MinimumBoundingGeometry_management(feature, "C:/Studia/PPG_II/Egzamin/Wyniki/rectangle_by_area.shp",
                                             "RECTANGLE_BY_AREA")
    arcpy.MinimumBoundingGeometry_management(feature, "C:/Studia/PPG_II/Egzamin/Wyniki/rectangle_by_width.shp",
                                             "RECTANGLE_BY_WIDTH")
    arcpy.MinimumBoundingGeometry_management(feature, "C:/Studia/PPG_II/Egzamin/Wyniki/convex_hull.shp", "CONVEX_HULL")
    arcpy.MinimumBoundingGeometry_management(feature, "C:/Studia/PPG_II/Egzamin/Wyniki/circle.shp", "CIRCLE")
    arcpy.MinimumBoundingGeometry_management(feature, "C:/Studia/PPG_II/Egzamin/Wyniki/envelope.shp", "ENVELOPE")

    lista = [r"C:/Studia/PPG_II/Egzamin/Wyniki/rectangle_by_area.shp",
             r"C:/Studia/PPG_II/Egzamin/Wyniki/rectangle_by_width.shp",
             r"C:/Studia/PPG_II/Egzamin/Wyniki/convex_hull.shp", r"C:/Studia/PPG_II/Egzamin/Wyniki/circle.shp",
             r"C:/Studia/PPG_II/Egzamin/Wyniki/envelope.shp"]

    lista_geometria_polygon = []
    lista_geometria_linie = []
    for i in lista:
        lista_geometria_linie.append(polygon_do_linii(i))
        for row in arcpy.da.SearchCursor(i, ["SHAPE@"]):
            lista_geometria_polygon.append(row[0])

    return lista_geometria_linie

def segment_length(punkt_1, punkt_2):
    x_1 = punkt_1.X
    y_1 = punkt_1.Y

    x_2 = punkt_2.X
    y_2 = punkt_2.Y

    dx = x_2 - x_1
    dy = y_2 - y_1

    return sqrt(dx ** 2 + dy ** 2)

def vertex_angle(punkt_1, punkt_2, punkt_3, obiekt, sprawdzenie):
    a = segment_length(punkt_1, punkt_2)
    b = segment_length(punkt_2, punkt_3)
    c = segment_length(punkt_1, punkt_3)

    if round(a+b, 10) == round(c, 10):
        kat = 180
    else:
        kat = acos((a ** 2 + b ** 2 - c ** 2) / (2 * a * b)) * 180 / pi

        x_1 = punkt_1.X
        y_1 = punkt_1.Y

        x_3 = punkt_3.X
        y_3 = punkt_3.Y

        punkt_srodek = arcpy.Point()
        punkt_srodek.X, punkt_srodek.Y = (x_1 + x_3) / 2, (y_1 + y_3) / 2

        if sprawdzenie == 'normalnie':
            if obiekt.disjoint(punkt_srodek):
                kat = 360 - kat
            else:
                pass
        else:
            if obiekt.disjoint(punkt_srodek):
                pass
            else:
                kat = 360 - kat

    return kat

def budynek(warstwa_wejsciowa, id_budynku):
    expression = "gmlID = " + "'" + id_budynku + "'"

    for row in arcpy.da.SearchCursor(warstwa_wejsciowa, ["OID@", "SHAPE@", "gmlId"], where_clause=expression):
        minimalne_geometrie = minimal_geometry(row[1])
        lista_czesci = []
        for part in row[1]:
            czesc = []
            for pnt in part:
                if pnt:
                    czesc.append(pnt)
                    # is_node(pnt, id_budynku, warstwa_wejsciowa)
                else:
                    trigger = True
                    lista_czesci.append(czesc)
                    czesc = []
            if len(lista_czesci) == 0 or trigger:
                trigger = False
                lista_czesci.append(czesc)

            tworzenie_warstwy_punktowej(lista_czesci)
            mdb = concave_hull("C:/Studia/PPG_II/Egzamin/Wyniki/Concave_Hull.shp", "C:/Studia/PPG_II/Egzamin/Wyniki/Concave_Hull_wynik.shp")
            minimalne_geometrie.append(mdb)
            wyniki = geometria_obiektu(lista_czesci, arcpy.Polygon(part), minimalne_geometrie)

    for i in wyniki:
        i.insert(0, id_budynku)

    return wyniki

def main():
    arcpy.env.overwriteOutput = 1

    warstwa_wejsciowa = r"C:/Studia/PPG_II/Egzamin/Dane/dane/BUBD.shp"

    id_budynkow = []
    for row in arcpy.da.SearchCursor(warstwa_wejsciowa, ["OID@", "SHAPE@", "gmlId"]):
        id_budynkow.append(row[2])

    wyniki = []
    for i in id_budynkow:
        wyniki.append(budynek(warstwa_wejsciowa, i))

    wynik_ostateczny = []
    for i in wyniki:
        for j in i:
            wynik_ostateczny.append(j)

    df = pd.DataFrame(np.array(wynik_ostateczny),
                      columns=['Identyfikator_budynku', 'Wierzcholek', 'Dlugosc_segmentu_przed', 'Dlugosc_segmentu_po',
                               'kat_wewnetrzny', 'Strzalka_do_boku_RECTANGLE_BY_AREA',
                               'Strzalka_do_boku_RECTANGLE_BY_WIDTH', 'Strzalka_do_boku_CONVEX_HULL',
                               'Strzalka_do_boku_CIRCLE', 'Strzalka_do_boku_ENVELOPE', 'Strzalka_do_boku_CONCAVE_HULL'])

    df.to_csv('C:/Users/Admin/Desktop/results.csv', index=False)

if __name__ == '__main__':
    main()
