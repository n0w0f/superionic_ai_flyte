from minio import Minio
from minio.error import S3Error
from dotenv import dotenv_values

# Load environment variables from .env file
env_vars = dotenv_values('../../.env')

# Retrieve MinIO endpoint, access key, and secret key from environment variables
endpoint = env_vars['MINIO_ENDPOINT']
access_key = env_vars['MINIO_ACCESS_KEY']
secret_key = env_vars['MINIO_SECRET_KEY']


# Set up MinIO client
minio_client = Minio(
    endpoint=endpoint,
    access_key=access_key,
    secret_key=secret_key,
    secure=False
)

print(minio_client.bucket_exists("superionic-ai"))



def upload_file(bucket_name, object_name, file_path):
    try:
        # Check if the bucket exists
        if not minio_client.bucket_exists(bucket_name):
            print(f"Bucket '{bucket_name}' does not exist.")
            return

        # Upload the file to MinIO
        minio_client.fput_object(bucket_name, object_name, file_path)

        print(f"File '{object_name}' uploaded successfully.")
    except S3Error as e:
        print(f"Error uploading file: {e}")

def download_file(bucket_name, object_name, file_path):
    try:
        # Download the file from MinIO
        minio_client.fget_object(bucket_name, object_name, file_path)

        print(f"File '{object_name}' downloaded successfully.")
    except S3Error as e:
        print(f"Error downloading file: {e}")

# Specify the bucket name, object name (file name in MinIO), and local file path
bucket_name = 'test-folder'
object_name = 'Na2PS3/mp-38200md_T600.log'
local_file_path = '/home/nawaf/workflows/superionic_ai_flyte/src/data/md_logs/Na2PS3/mp-38200md_T600.log'

# Upload the file to MinIO
upload_file(bucket_name, object_name, local_file_path)

# Download the file from MinIO
#download_file(bucket_name, object_name = 'Na2PS3/mp-38200md_T600.log',file_path ="/home/nawaf/workflows/superionic_ai_flyte/src/data/md_logs/Na2PS3")
