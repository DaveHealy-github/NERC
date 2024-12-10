classdef TriTransferFunction
   
    % DESCRIPTION:
    %   Base Transfer function class that supports multi-phase flow
    
    properties
        nphases
    end
    methods
        function transferfunction = TriTransferFunction()
            
        end
        
        function [] = validate_fracture_matrix_structures(tf,fracture_fields,matrix_fields, matrix2_fields)

            fracture_fields = fieldnames(fracture_fields);
            matrix_fields = fieldnames(matrix_fields);
            matrix2_fields = fieldnames(matrix2_fields);
            
            if(tf.nphases == 1)
                assert(any(ismember(fracture_fields,'pf')),...
                      'Could not find "pof" field in "fracture_fields" structure')
                  
                assert(any(ismember(matrix_fields,'pm')),...
                    'Could not find "pom" field in "matrix_fields" structure')

                assert(any(ismember(matrix2_fields,'pm2')),...
                    'Could not find "pom2" field in "matrix2_fields" structure')
            end
            if(tf.nphases == 2)
                assert(any(ismember(fracture_fields,'pof')),...
                    'Could not find "pof" field in "fracture_fields" structure')
                assert(any(ismember(fracture_fields,'swf')),...
                    'Could not find "swf" field in "fracture_fields" structure')
                
                assert(any(ismember(matrix_fields,'pom')),...
                    'Could not find "pom" field in "matrix_fields" structure')
                assert(any(ismember(matrix_fields,'swm')),...
                    'Could not find "swm" field in "matrix_fields" structure')

                assert(any(ismember(matrix2_fields,'pom2')),...
                    'Could not find "pom2" field in "matrix2_fields" structure')
                assert(any(ismember(matrix2_fields,'swm2')),...
                    'Could not find "swm2" field in "matrix2_fields" structure')
            end
            if(tf.nphases == 3) % 3 phases
                assert(any(ismember(fracture_fields,'pof')),...
                    'Could not find "pof" field in "fracture_fields" structure')
                assert(any(ismember(fracture_fields,'swf')),...
                    'Could not find "swf" field in "fracture_fields" structure')
                assert(any(ismember(fracture_fields,'sgf')),...
                    'Could not find "sgf" field in "fracture_fields" structure')
                
                assert(any(ismember(matrix_fields,'pom')),...
                    'Could not find "pom" field in "matrix_fields" structure')
                assert(any(ismember(matrix_fields,'swm')),...
                    'Could not find "swm" field in "matrix_fields" structure')
                assert(any(ismember(matrix_fields,'sgm')),...
                    'Could not find "sgf" field in "matrix_fields" structure')

                assert(any(ismember(matrix2_fields,'pom2')),...
                    'Could not find "pom2" field in "matrix2_fields" structure')
                assert(any(ismember(matrix2_fields,'swm2')),...
                    'Could not find "swm2" field in "matrix2_fields" structure')
                assert(any(ismember(matrix2_fields,'sgm2')),...
                    'Could not find "sgf2" field in "matrix2_fields" structure')
            end
            
        end
        
        function [Gamma] = calculate_transfer(tf,model,fracture_fields,matrix_fields, matrix2_fields)
            % This method should be reimplemented in the derived classes
            Gamma{1} = 0;
            Gamma{2} = 0;
            Gamma{3} = 0;
        end
        
    end
    
    
end
