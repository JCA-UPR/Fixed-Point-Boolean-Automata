<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <link rel="stylesheet" href="style.css">
	<title>File Upload</title>
</head>
<body>
    <h2>Automata Website</h2>
    <form action="upload.php" method="post" enctype="multipart/form-data">
    <input type="file" name="file">
    <button type="submit" name="submit">upload</button>
    </form>
    
    <a href="./output.txt" download> <button>Download File</button>
    </a>
    <?php
	$result = shell_exec("python3 main.py ./uploads/input.csv");
	echo "PHP got the result $result";
	?>
</body>
</html>
