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
library(ggplot2)

# Function Definitions ----------------------------------------------------

extract_transform_load <- function(filename, threshold = 5000) {
  my_tbl <- readr::read_csv(filename) %>%
    dplyr::select(
      well = `Well Name`,
      count = `All Events #Events`,
      mean_fluorescence = `All Events FITC-A Mean`
    ) %>%
    dplyr::mutate(
      plate = factor(rep(c("A", "B"), times = 3, each = 96)),
      ligand = factor(rep(c("None", "VAC", "Vanillin"), times = 2, each = 96)),
      type = "Mutant",
      type = factor(replace(type, well %in% c("H11", "H12"), "Control")),
    )
    # ) %>%
    # dplyr::filter(count >= threshold)
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
  scale_fill_brewer(palette = "Set1")

ggplot(my_df, aes(x = type,
                  y = mean_fluorescence,
                  color = type)) +
  geom_boxplot() +
  facet_grid(. ~ ligand) +
  scale_color_brewer(palette = "Set1", guide = FALSE)

reference <- filter(my_df, ligand == "None")$mean_fluorescence

my_df <- my_df %>%
  dplyr::group_by(ligand) %>%
  dplyr::mutate(
    fc = mean_fluorescence / reference
  )

ggplot(my_df, aes(x = fc,
                  fill = type)) +
  geom_histogram(position = "dodge") +
  facet_grid(. ~ ligand) +
  scale_fill_brewer(palette = "Set1")
