from minio import Minio
from minio.error import S3Error
from dotenv import dotenv_values

from typing import List
import os 
import io

# Load environment variables from .env file
env_vars = dotenv_values('/config/.env')

# Retrieve MinIO endpoint, access key, and secret key from environment variables
endpoint = '192.168.9.13:30002'
access_key = 'GbhgmzTziWIanCKPaay0'
secret_key = 'miniostorage'

# Use the variables as needed
print(f"MinIO Endpoint: {endpoint}")


# Set up MinIO client
remote_client = Minio(
    endpoint=endpoint,
    access_key=access_key,
    secret_key=secret_key,
    secure=False
)

print(remote_client.bucket_exists("superionic-ai"))

def upload_file( file_path : str, config_local : str, config_remote : str , bucket_name :str ):

    try:
        # Check if the bucket exists
        if not remote_client.bucket_exists(bucket_name):
            print(f"Bucket '{bucket_name}' does not exist.")
            return
 
        # Generate the MinIO object name by appending the file name to a base folder path
        object_name =  file_path.replace( config_local, config_remote )


        # Upload the file to MinIO

        remote_client.fput_object(bucket_name, object_name, file_path)


        print(f"File '{file_path}' uploaded as '{object_name}'")

        return object_name
        
    except S3Error as e:
        print(f"Error uploading files to MinIO: {e}")




#object name is minio path
def download_file(object_name : str, file_path :str, bucket_name = "superionic-ai"):
    try:
        # Download the file from MinIO
        remote_client.fget_object(bucket_name, object_name, file_path)

        print(f"File '{object_name}' downloaded successfully.")
    except S3Error as e:
        print(f"Error downloading file: {e}")




def upload_files_to_remote( file_paths :str, config_local :str, config_remote :str , bucket_name  = "superionic-ai"):
    try:
        # Check if the bucket exists
        if not remote_client.bucket_exists(bucket_name):
            print(f"Bucket '{bucket_name}' does not exist.")
            return

        for file_path in file_paths:
            
            # Generate the MinIO object name by appending the file name to a base folder path
            object_name =  file_path.replace( config_local, config_remote )


            # Upload the file to MinIO

            remote_client.fput_object(bucket_name, object_name, file_path)


            print(f"File '{file_path}' uploaded as '{object_name}'")
        
        print("All files uploaded successfully.")
    except S3Error as e:
        print(f"Error uploading files to MinIO: {e}")



def download_folder_from_remote( folder_name : str  ,   destination_path :str,  bucket_name,):
    try:
        # Check if the bucket exists
        if not remote_client.bucket_exists(bucket_name):
            print(f"Bucket '{bucket_name}' does not exist.")
            return

        # Create the destination folder if it doesn't exist
        os.makedirs(destination_path, exist_ok=True)

        # Retrieve a list of objects in the folder
        objects = remote_client.list_objects(bucket_name, prefix=folder_name, recursive=True)


        for obj in objects:

            # Skip directories
            if obj.is_dir:
                continue

            # Generate the local file path by replacing the MinIO folder name with the destination path
   
            local_file_path = obj.object_name.replace(folder_name, destination_path)

            # Create the directory structure if it doesn't exist
            os.makedirs(os.path.dirname(local_file_path), exist_ok=True)

            # Download the file from MinIO
            remote_client.fget_object(bucket_name, obj.object_name, local_file_path)
     


            print(f"File '{obj.object_name}' downloaded to '{local_file_path}'")

        print("Folder downloaded successfully.")
    except S3Error as e:
        print(f"Error downloading folder from MinIO: {e}")




def create_folders_remote (bucket_name : str , elements: List[str], root_path: str):

    try:
        # Check if the bucket exists
        if not remote_client.bucket_exists(bucket_name):
            print(f"Bucket '{bucket_name}' does not exist.")
            return

        # Create an empty object with a key representing the folder
        for element in elements:

            folder_path = os.path.join(root_path, element)
            object_name = folder_path + '/'
            data = b''
            data_stream = io.BytesIO(data)
            remote_client.put_object(bucket_name, object_name, data_stream, 0)

            print(f"Folder '{object_name}' created successfully.")
    except S3Error as e:
        print(f"Error creating folder: {e}")



