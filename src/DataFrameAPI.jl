### -*- Mode: Julia -*-

"Module DataFrameAPI

"
module DataFrameAPI

using Graphs
using CSV
using MetaGraphs
using DataFrames

function find_edges_from_vertex(path_dataframe_edges::String, first_vertex::Int64, second_vertex::Int64)
    df_edges_file = CSV.File(path_dataframe_edges)
    df_edges = DataFrame(df_edges_file)
    # Apply filter
    filtered_edges = filter(row -> (row.origin == first_vertex && row.destination == second_vertex) ||
                                   (row.origin == second_vertex && row.destination == first_vertex), df_edges)
    new_df_edges = DataFrame(filtered_edges)
    return new_df_edges
end

function find_neighbors_from_vertex(path_dataframe_edges::String, path_dataframe_labels::String, vertex::Int64)
    df_edges_file = CSV.File(path_dataframe_edges)
    df_edges = DataFrame(df_edges_file)
    df_labels_file = CSV.File(path_dataframe_labels)
    df_labels = DataFrame(df_labels_file)
    # Extract neighbors vertex based on df_edges
    neighbors = unique(vcat(filter(row -> row.origin == vertex, df_edges).destination,
                            filter(row -> row.destination == vertex, df_edges).origin))
    # Apply filter row in df_labels based on label field
    neighbor_labels = filter(row -> row.label in neighbors, df_labels)
    return neighbor_labels
end

function find_vertex_from_CI(path_dataframe_labels::String, CI::CartesianIndex)
    df_labels_file = CSV.File(path_dataframe_labels)
    df_labels = DataFrame(df_labels_file)
    position = string(CI)
    filtered_labels = filter(row -> row.position_label == position, df_labels)
    df_label_result = DataFrame(filtered_labels)
    return df_label_result
end

function find_vertex_from_color(path_dataframe_labels::String, color::Float64)
    df_labels_file = CSV.File(path_dataframe_labels)
    df_labels = DataFrame(df_labels_file)
    filtered_labels = filter(row -> row.color_label == color, df_labels)
    df_label_result = DataFrame(filtered_labels)
    return df_label_result
end

function find_vertex_from_label(path_dataframe_labels::String, label::Int64)
    df_labels_file = CSV.File(path_dataframe_labels)
    df_labels = DataFrame(df_labels_file)
    filtered_labels = filter(row -> row.label == label, df_labels)
    df_label_result = DataFrame(filtered_labels)
    return df_label_result
end

function count_vertex_degree(path_dataframe_labels::String, vertex::Int64)
    df_labels_file = CSV.File(path_dataframe_labels)
    df_labels = DataFrame(df_labels_file)
    outgoing_degree = nrow(filter(row -> row.origin == vertex, df_labels))
    incoming_degree = nrow(filter(row -> row.destination == vertex, df_labels))
    return outgoing_degree + incoming_degree
end

function count_edges(path_dataframe_edges::String)
    df_edges_file = CSV.File(path_dataframe_edges)
    df_edges = DataFrame(df_edges_file)
    num_records = nrow(df_edges)
    return num_records
end

function count_vertices(path_dataframe_labels::String)
    df_labels_file = CSV.File(path_dataframe_labels)
    df_labels = DataFrame(df_labels_file)
    num_records = nrow(df_labels)
    return num_records
end

end # module

### end of file -- DataFrameAPI.jl
