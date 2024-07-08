

configfile: "/lustre/pulkit.h/snakemake_local/snakemake_tutorial/configs/config.yaml" # configuration file



rule all:
    input:
        sample_out = f"{config['home_dir']}/output/RDS_File/RDS_File.rds",
        filtering_out = f"{config['home_dir']}/output/RDS_File/Filtered_data.rds",
        clustering_out1 =f"{config['home_dir']}/output/csv/Meta_Data.csv",
        clustering_out2 = f"{config['home_dir']}/output/RDS_File/Clustering.rds",
        clustering_file3 = f"{config['home_dir']}/output/csv/Clustering_Resolutions.csv",
        find_markers_out1 = f"{config['home_dir']}/output/csv/All_markers.txt",
        find_markers_out2 = f"{config['home_dir']}/output/csv/Top_Markers.txt",

rule sample:
    input:
        file= f"{config['home_dir']}/input/matrix/SC3pv3_GEX_Human_PBMC_filtered_feature_bc_matrix.h5",
    output:
        file= f"{config['home_dir']}/output/RDS_File/RDS_File.rds",
        
    params:
        file_no = 3,
        no_features= 200,
        pdf1 = f"{config['home_dir']}/output/pdfs/Number_of_Genes.pdf",
        pdf2 = f"{config['home_dir']}/output/pdfs/Miochondrial_percentage.pdf",
        pdf3 = f"{config['home_dir']}/output/pdfs/scatter_plot.pdf",
        
    shell:
        """
        echo "Filter out features that are present in fewer than 'n' cells (n)=  "
        read no_of_cells
        echo "Filter out cells which have features less than < "
        read no_of_features
        Rscript scripts/sample.R {input.file} {output.file} $no_of_cells $no_of_features {params.pdf1} {params.pdf2} {params.pdf3}
        """


rule filtering:
    input:
        file= f"{config['home_dir']}/output/RDS_File/RDS_File.rds",
    output:
        file= f"{config['home_dir']}/output/RDS_File/Filtered_data.rds", 
    params:
        # no_of_features_grt = 100,
        # no_of_features_less = 6000, 
        # mitochondrial_percent = 25,
        # find_variables = 2000,
        # top_features = 10,
        pdf2 = f"{config['home_dir']}/output/pdfs/Miochondrial_percentage.pdf",
        pdf4 = f"{config['home_dir']}/output/pdfs/Highly_variable_Genes.pdf",
        pdf5 = f"{config['home_dir']}/output/pdfs/No_of_Cells_and_Features.pdf",
        pdf6 = f"{config['home_dir']}/output/pdfs/Scaled_Highly_variable_Genes.pdf",
        pdf7 = f"{config['home_dir']}/output/pdfs/PCA.pdf",
        pdf8 = f"{config['home_dir']}/output/pdfs/Top_10_PC.pdf" ,
        pdf9 = f"{config['home_dir']}/output/pdfs/Elbow_plot.pdf" ,
        
    shell:
        """
        code {params.pdf2}
        echo "Select minium number of Features count = "
        read no_of_features_grt
        echo "Select Maximum number of Features count  = "
        read no_of_features_less
        echo "Filter out Cells which have Mitochondrial Percentage Less Than = "
        read mitochondrial_percent
        echo "Number of Features which are Highly Variable = "
        read variable_features
        echo "Number of Features to be Highlighted = "
        read top_features
        code {params.pdf6}
        Rscript scripts/Filtering_sc_data.R {input.file} {output.file} $no_of_features_grt $no_of_features_less $mitochondrial_percent $variable_features $top_features {params.pdf4} {params.pdf5} {params.pdf6} {params.pdf7} {params.pdf8} {params.pdf9} 
        """


rule clustering:
    input:
        file = f"{config['home_dir']}/output/RDS_File/Filtered_data.rds"
    output:
        file  = f"{config['home_dir']}/output/csv/Meta_Data.csv",
        file2 = f"{config['home_dir']}/outputRDS_File/Clustering.rds",
        file3 = f"{config['home_dir']}/output/csv/Clustering_Resolutions.csv"
    params:
        pdf9 =f"{config['home_dir']}/output/pdfs/Elbow_plot.pdf" ,
    shell:
        """
        code {params.pdf9}
        echo "Please enter the number of PCS dimensions to use for finding clusters: "
        read dims
        Rscript scripts/Clustering_of_SC.R {input.file} {output.file} {output.file2} $dims {output.file3}
        code {output.file3}
        """




rule find_markers:
    input:
        file= f"{config['home_dir']}/output/RDS_File/Clustering.rds"
    output:
        file = f"{config['home_dir']}/output/csv/All_markers.txt",
        file2= f"{config['home_dir']}/output/csv/Top_Markers.txt"
    params:
        Onlypositive_FC_= "TRUE",
        Resolutions= "RNA_snn_res.0.2",
        min_pct= 0.25,
        Top_listed_Genes= 10,
        pdf1= f"{config['home_dir']}/output/pdfs/Heatmap_of_markers.pdf",
        pdf2= f"{config['home_dir']}/output/pdfs/Clusters_of_Cells.pdf",
        dims= 15  # Default value for dims
    shell:
        """
        echo "Please enter the number of Resolutions: "
        read Resolutions
        code {params.pdf1}
        Rscript scripts/Find_Markers.R {input.file} {output.file} {output.file2} $Resolutions {params.Onlypositive_FC_} {params.min_pct} {params.Top_listed_Genes} {params.pdf1} {params.pdf2} {params.dims}
        code {params.pdf2}
        """

