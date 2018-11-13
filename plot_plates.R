# Copyright 2018 Novo Nordisk Foundation Center for Biosustainability,
# Technical University of Denmark.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# Libraries ---------------------------------------------------------------

library(readr)
library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)

# Function Definitions ----------------------------------------------------

extract_transform_load <- function(filename, threshold = 5000) {
  # Load the tabular data file.
  my_tbl <- readr::read_csv(filename) %>%
    # Select specific columns and rename them to be simpler.
    dplyr::select(
      well = `Well Name`,
      count = `All Events #Events`,
      mean_fluorescence = `All Events FITC-A Mean`,
      time = `Record Date`
    ) %>%
    # Create additional columns based on the a priori knowledge about the data
    # layout.
    dplyr::mutate(
      time = lubridate::mdy_hms(time),
      plate = factor(rep(c("A", "B"), times = 3, each = 96)),
      ligand = factor(rep(c("None", "VAC", "Vanillin"), each = 192)),
      type = "Mutant",
      type = factor(replace(type, well %in% c("H11", "H12"), "Control")),
      position = interaction(plate, well)
    )
  # Compute fold changes as compared to the reference without any ligand.
  fold_change <- dplyr::select(my_tbl, position, ligand, mean_fluorescence) %>%
    tidyr::spread(ligand, mean_fluorescence) %>%
    dplyr::mutate(
      VAC = VAC / None,
      Vanillin = Vanillin / None,
      None = None / None
    ) %>%
    tidyr::gather(None, VAC, Vanillin, key = "ligand", value = "fc")
  # Join the two tables again based on common position and ligand.
  my_tbl <- dplyr::inner_join(my_tbl, fold_change) %>%
    # Remove all rows with too low count.
    dplyr::filter(count >= threshold)

  return(my_tbl)
}

# Command Line Parsing ----------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Expected one command line argument specifying the CSV data file.")
}

# Load Data ---------------------------------------------------------------

my_df <- extract_transform_load(args[1])

# Plot Data ---------------------------------------------------------------

ggplot(my_df, aes(x = mean_fluorescence,
                  fill = type)) +
  geom_histogram(position = "dodge") +
  facet_grid(. ~ ligand) +
  scale_fill_brewer(palette = "Set1") +
  coord_cartesian(xlim = c(0, 8000))

ggplot(my_df, aes(x = type,
                  y = mean_fluorescence,
                  color = type)) +
  geom_boxplot() +
  facet_grid(. ~ ligand) +
  scale_color_brewer(palette = "Set1", guide = FALSE) +
  coord_cartesian(ylim = c(0, 8000))

ggplot(filter(my_df, ligand != "None"), aes(x = fc, fill = type)) +
  geom_histogram(position = "dodge") +
  facet_grid(. ~ ligand) +
  scale_fill_brewer(palette = "Set1")

# Low VAC fold changes.
low_vac <- dplyr::filter(my_df, ligand == "VAC", fc < 1.15)
vanillin <- dplyr::filter(my_df, ligand == "Vanillin", position %in% low_vac$position)
controls <- dplyr::filter(my_df, type == "Control", ligand != "None")
variants <- dplyr::bind_rows(low_vac, vanillin, controls)

ggplot(variants, aes(x = time, y = fc, color = type)) +
  geom_point() +
  facet_grid(. ~ ligand) +
  scale_color_brewer(palette = "Set1")
