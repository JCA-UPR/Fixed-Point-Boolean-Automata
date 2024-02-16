<?php
echo shell_exec("python script.py");
if (isset($_POST['submit'])) {
    $file = $_FILES['file'];

    $fileName = $file['name'];
    $fileTmpName = $file['tmp_name'];
    $fileSize = $file['size'];
    $fileError = $file['error'];

    // Check if there is no error
    if ($fileError === 0) {
        $uploadPath = 'uploads/input.csv'; // Directory where to upload the file
        echo $fileName;
        move_uploaded_file($fileTmpName, $uploadPath);
        echo "File uploaded successfully!";
        echo $_FILE['file'];
    } else {
        echo "Error uploading file.";
    }
}

?>
<!DOCTYPE html>
<html>
<head>
    <title>Go Back</title>
</head>
<body>
    <h2>Click the button to go back</h2>
    <form method="post">
        <button type="submit" name="goback">Go Back</button>
    </form>
</body>
</html>

<?php
if (isset($_POST['goback'])) {
    header('Location: index.php');
    exit();
}
?>
