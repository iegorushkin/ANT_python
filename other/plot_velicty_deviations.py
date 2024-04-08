from functions_visualization import (extract_coastlines,
                                     visualize_velocity_deviations,
                                     visualize_velocity_deviations_ru)

# paths
coastlines_path = "E:\\Work\\Projects\\Central_Kamchatka\\Temp\\coastal_line_4.bln"
stations_path = "E:\\Work\\Projects\\Central_Kamchatka\\Temp\\stations_locations.xlsx"
rays_path = "E:\\Work\\Projects\\Central_Kamchatka\\Temp\\rays"
save_path="E:\\Work\\Projects\\Central_Kamchatka\\_Paper\\Figures"

## Coordinates for coastlines extraction
# longitude
x_min = 155.6
x_max = 160.4
# latitude
y_min = 52.4
y_max = 54.6

## Coordinates of the boundaries
## of the group velocity deviation visualization maps
# longitude
x_limits = [155.8, 160.2]
# latitude
y_limits = [52.6, 54.45]

coastlines_list = extract_coastlines(coastlines_path, x_min, x_max, y_min, y_max)
visualize_velocity_deviations(coastlines_list, x_limits, y_limits,
                              stations_path, rays_path,
                              step=3, yx_aspect_ratio=1.68, key=1,
                              dont_open=1, save_path=save_path)
# visualize_velocity_deviations_ru(coastlines_list, x_limits, y_limits,
#                                  stations_path, rays_path,
#                                  step=3, yx_aspect_ratio=1.68, key=1,
#                                  dont_open=1, save_path=save_path, ru=True)
