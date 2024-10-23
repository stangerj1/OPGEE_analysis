# %%
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point, Polygon, MultiPolygon
from fuzzywuzzy import process, fuzz
from geopy.distance import geodesic

# Load GEM data
data = pd.read_csv("GEM_new.csv", encoding='latin1')
columns_to_drop = ['Subnational unit (province, state)', 'Status year', 'Discovery year', 
                   'Production start year', 'Basin', 'Concession / block']
data = data.drop(columns=columns_to_drop)
data['normalized_field_name'] = data['Unit name'].str.lower().str.strip()

# Load Woodmac data
df = pd.read_csv("Woodmac_for_match.csv", encoding='latin1')
df = df[df['field_is_top_level'] == 'Y']
df['normalized_field_name'] = df['field_name'].str.lower().str.strip()

# Create GeoDataFrame for GEM data
gem_gdf = gpd.GeoDataFrame(data, geometry=gpd.points_from_xy(data.Longitude, data.Latitude))
gem_gdf.set_crs(epsg=4326, inplace=True)

# Load Woodmac polygon data
gdf = gpd.read_file("field_geospatial_summary.gpkg")
columns_to_extract_WM = ['field_name','field_group','field_is_discovery','country_name','field_is_top_level',
                         'region','onshore_offshore_tags','basin_name','sector_name','sector_country',
                         'field_onshore_offshore','field_status','field_year_discovery','field_year_production_start','field_operator',
                         'field_is_parent','field_parent','field_taxation','field_centroid_x','field_centroid_y','geometry']
woodmac_polygons = gdf[columns_to_extract_WM]
woodmac_polygons.set_crs(epsg=4326, inplace=True)

# %%
# Perform sublevel polygon matching
sublevel_polygons = woodmac_polygons[woodmac_polygons['field_is_top_level'] == 'N']
sublevel_polygon_matches = []

for idx, gem_row in gem_gdf.iterrows():
    for poly_idx, woodmac_row in sublevel_polygons.iterrows():
        if woodmac_row.geometry.contains(gem_row.geometry):
            combined_row = {'field_name': woodmac_row['field_name'],
                            'field_group': woodmac_row['field_group'],
                            'field_centroid_y': woodmac_row['field_centroid_y'],
                            'field_centroid_x': woodmac_row['field_centroid_x'],
                             **gem_row.to_dict(),
                            'polygon_geometry': woodmac_row.geometry,}
            sublevel_polygon_matches.append(combined_row)

sublevel_polygon_matched_fields = pd.DataFrame(sublevel_polygon_matches)

# %%
print(f"Number of matched fields: {len(sublevel_polygon_matched_fields)}")
matched_unit_names = sublevel_polygon_matched_fields['Unit name'].unique()
unmatched_gem_gdf = gem_gdf[~gem_gdf['Unit name'].isin(matched_unit_names)]

# %%
toplevel_polygons = woodmac_polygons[woodmac_polygons['field_is_top_level'] == 'Y']
toplevel_polygon_matches = []

for idx, gem_row in unmatched_gem_gdf.iterrows():
    for poly_idx, woodmac_row in toplevel_polygons.iterrows():
        if woodmac_row.geometry.contains(gem_row.geometry):
            combined_row = {'field_name': woodmac_row['field_name'],
                            'field_group': woodmac_row['field_group'],
                            'field_centroid_y': woodmac_row['field_centroid_y'],
                            'field_centroid_x': woodmac_row['field_centroid_x'],
                             **gem_row.to_dict(),
                            'polygon_geometry': woodmac_row.geometry,}
            toplevel_polygon_matches.append(combined_row)

toplevel_polygon_matched_fields = pd.DataFrame(toplevel_polygon_matches)
polygon_matched_fields= []
polygon_matched_fields = pd.concat([sublevel_polygon_matched_fields, toplevel_polygon_matched_fields]).drop_duplicates(subset=['Unit name']).reset_index(drop=True)

matched_unit_names = polygon_matched_fields['Unit name'].unique()
unmatched_data = gem_gdf[~gem_gdf['Unit name'].isin(matched_unit_names)]

# %%
# Exact name match
exact_matches = []
for idx, row in unmatched_data.iterrows():
    exact_match_row = df[df['normalized_field_name'] == row['normalized_field_name']]
    if not exact_match_row.empty:
        exact_match_row = exact_match_row.iloc[0]
        combined_row = {'field_name': exact_match_row['field_name'],
                        'field_group': exact_match_row['field_group'],
                        'field_centroid_y': exact_match_row['field_centroid_y'],
                        'field_centroid_x': exact_match_row['field_centroid_x'],
                        **row.to_dict()}
        exact_matches.append(combined_row)

exact_matched_fields = pd.DataFrame(exact_matches)

unmatched_data = unmatched_data[~unmatched_data['Unit name'].isin(exact_matched_fields['Unit name'])]

# %%
# Fuzzywuzzy match
def get_best_fuzzy_match(row, choices, threshold=89):
    match, score = process.extractOne(row['normalized_field_name'], choices)
    return match if score >= threshold else None

choices = df['normalized_field_name'].tolist()

unmatched_data['fuzzy_match'] = unmatched_data.apply(
    lambda row: get_best_fuzzy_match(row, choices), axis=1
)
df['matched'] = False
fuzzy_matches = []
fuzzy_matches_without_location = []

for idx, row in unmatched_data.iterrows():
    if row['fuzzy_match'] is not None:
        matched_row = df[(
            df['normalized_field_name'] == row['fuzzy_match']
            ) & (df['matched'] == False)]
        if not matched_row.empty:
            matched_row = matched_row.iloc[0]
            if pd.isna(row['Latitude']) or pd.isna(row['Longitude']):
                combined_row = {'field_name': matched_row['field_name'],
                                'field_group': matched_row['field_group'],
                                'field_centroid_y': matched_row['field_centroid_y'],
                                'field_centroid_x': matched_row['field_centroid_x'],
                                **row.to_dict()}
                fuzzy_matches_without_location.append(combined_row)
                df.loc[matched_row.name, 'matched'] = True
            else:
                distance = geodesic((row['Latitude'], row['Longitude']),
                                    (matched_row['field_centroid_y'], matched_row['field_centroid_x'])).kilometers
                if distance <= 20:
                    combined_row = {'field_name': matched_row['field_name'],
                                    'field_group': matched_row['field_group'],
                                    'field_centroid_y': matched_row['field_centroid_y'],
                                    'field_centroid_x': matched_row['field_centroid_x'],
                                     **row.to_dict()}
                    fuzzy_matches.append(combined_row)
                    df.loc[matched_row.name, 'matched'] = True

fuzzy_matched_fields = pd.DataFrame(fuzzy_matches)
fuzzy_matches_without_location_df = pd.DataFrame(fuzzy_matches_without_location)

unmatched_data = unmatched_data[~unmatched_data['Unit name'].isin(fuzzy_matched_fields['Unit name'])]

# %%
# Spatial match
def filter_fuzzy_match(row, target_row, threshold=41):
    score = fuzz.ratio(row['normalized_field_name'], target_row['normalized_field_name'])
    return score >= threshold

def find_closest_match(row, df, max_distance=8.0, threshold=41):
    closest_match = None
    closest_distance = max_distance
    if pd.isna(row['Latitude']) or pd.isna(row['Longitude']):
        return closest_match
    row_location = (row['Latitude'], row['Longitude'])
    
    for idx, potential_match in df.iterrows():
        if potential_match['matched']:
            continue
        if pd.isna(potential_match['field_centroid_y']) or pd.isna(potential_match['field_centroid_x']):
            continue
        potential_location = (potential_match['field_centroid_y'], potential_match['field_centroid_x'])
        distance = geodesic(row_location, potential_location).kilometers
        if distance < closest_distance and filter_fuzzy_match(row, potential_match, threshold):
            closest_distance = distance
            closest_match = potential_match
    return closest_match

df['matched'] = False

matched_records = []
for idx, row in unmatched_data.dropna(subset=['Latitude', 'Longitude']).iterrows():
    closest_match = find_closest_match(row, df, max_distance=8.0)
    if closest_match is not None:
        combined_row = {'field_name': closest_match['field_name'], 
                        'field_group': closest_match['field_group'],
                        'field_centroid_y': closest_match['field_centroid_y'],
                        'field_centroid_x': closest_match['field_centroid_x'],
                        **row.to_dict()}
        matched_records.append(combined_row)
        df.loc[closest_match.name, 'matched'] = True

spatial_matched_fields = pd.DataFrame(matched_records)

# %%
all_matches = pd.concat([polygon_matched_fields, exact_matched_fields, fuzzy_matched_fields, fuzzy_matches_without_location_df, spatial_matched_fields]).drop_duplicates(subset=['Unit name']).reset_index(drop=True)
print(f"Number of matched fields: {len(all_matches)}")
print(all_matches[['Unit name', 'field_name', 'Latitude', 'Longitude', 'field_centroid_y', 'field_centroid_x', 'field_group']].head())

all_matches.to_csv("global_polygon_all.csv", index=False)

# %%
print(f"Number of matched fields: {len(fuzzy_matched_fields)}")

# %%
print(f"Number of matched fields: {len(spatial_matched_fields)}")

# %%
# find nearest polygon for unmatched oil field
unmatched_data = data[~data['Unit name'].isin(all_matches['Unit name'])]
woodmac_polygons['normalized_field_name'] = woodmac_polygons['field_name'].str.lower().str.strip()
def filter_fuzzy_match_30(row, target_row, threshold=30):
    score = fuzz.ratio(row['normalized_field_name'], target_row['normalized_field_name'])
    return score >= threshold

nearest_polygon_matches = []

for idx, row in unmatched_data.iterrows():
    if pd.isna(row['Latitude']) or pd.isna(row['Longitude']):
        continue
    row_location = (row['Latitude'], row['Longitude'])
    closest_polygon = None
    closest_distance = float('inf')
    
    for poly_idx, woodmac_row in woodmac_polygons.iterrows():
        polygon = woodmac_row.geometry
        centroid = polygon.centroid
        distance = geodesic(row_location, (centroid.y, centroid.x)).kilometers
        if distance < closest_distance and filter_fuzzy_match_30(row, woodmac_row, threshold=30):
            closest_distance = distance
            closest_polygon = woodmac_row
    
    if closest_polygon is not None:
        combined_row = {'field_name': closest_polygon['field_name'],
                        'field_group': closest_polygon['field_group'],
                        'field_centroid_y': closest_polygon['field_centroid_y'],
                        'field_centroid_x': closest_polygon['field_centroid_x'],
                        **row.to_dict(),
                        'polygon_geometry': closest_polygon.geometry,
                        'distance_to_polygon': closest_distance}
        nearest_polygon_matches.append(combined_row)

nearest_polygon_matched_fields = pd.DataFrame(nearest_polygon_matches)

# %%
nearest_polygon_matched_fields = nearest_polygon_matched_fields.sort_values(by='distance_to_polygon').reset_index(drop=True)

print(f"Number of nearest polygon matched fields: {len(nearest_polygon_matched_fields)}")
print(nearest_polygon_matched_fields[['Unit name', 'field_name', 'Latitude', 'Longitude', 'field_centroid_y', 'field_centroid_x', 'field_group', 'distance_to_polygon']].head())

nearest_polygon_matched_fields.to_csv("global_nearest_polygon_all.csv", index=False)

# %%



